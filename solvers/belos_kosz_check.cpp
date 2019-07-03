#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosBlockCGSolMgr.hpp>

#include <Tpetra_MatrixIO.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Operator.hpp>
#include <Kokkos_DefaultNode.hpp>

#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_BorderedOperator.hpp>

#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>

#include "../evolution/matrices.h"
#include <random>
#include <cstring>
#include <sys/stat.h>
#include <experimental/filesystem>

#include "/home/kdeweese/TreePCG/haoran_code/io.h"
#include "/home/kdeweese/TreePCG/haoran_code/graph.h"
#include "/home/kdeweese/TreePCG/haoran_code/kosz.h"

using namespace Teuchos;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using std::endl;
using std::cout;
using std::vector;
using Teuchos::tuple;


typedef double                        ST;
typedef ScalarTraits<ST>                SCT;
typedef SCT::magnitudeType               MT;
typedef Tpetra::Operator<ST,int>         OP;
typedef Tpetra::MultiVector<ST,int>      MV;
typedef Belos::OperatorTraits<ST,MV,OP> OPT;
typedef Belos::MultiVecTraits<ST,MV>    MVT;
typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
typedef Tpetra::Map<int, int, Node>     Map;
typedef Tpetra::CrsMatrix<ST,int,int> TCrsMatrix;

struct result{
  ST residual;
  double time;
  int iters;
  int success;
  int flops;
  std::string filename;
  bool operator>(const result& rhs) const {
    if(flops == rhs.flops)
      return residual > rhs.residual;
    else
      return flops > rhs.flops;
  }
};


class SingularOp : public virtual OP {
public:
  // Constructor
  SingularOp (RCP<TCrsMatrix> matrix, RCP<MV> nullspace) : matrix_(matrix),nullspace_(nullspace) {}

  // Destructor
  virtual ~SingularOp() {}

  RCP<const Map> getDomainMap() const { return matrix_->getDomainMap(); }
  RCP<const Map> getRangeMap() const { return matrix_->getRangeMap(); }

  // Computes Y = alpha Op X + beta Y
  void apply (const MV& X, MV& Y,
              Teuchos::ETransp mode = Teuchos::NO_TRANS,
              ST alpha = SCT::one (),
              ST beta = SCT::zero ()) const
  {
    matrix_->apply(X,Y,mode,alpha,beta);
    

    std::vector<ST> product1(1);
    std::vector<ST> product2(1);
    MVT::MvDot(*nullspace_, Y, product1);
    MVT::MvDot(*nullspace_, *nullspace_, product2);
    MVT::MvAddMv(-product1[0]/product2[0],*nullspace_,1.0,Y,Y);
  
    
  }

  // Whether the operator supports applying the transpose
  bool hasTransposeApply() const { return true; }

private:
  RCP<TCrsMatrix> matrix_;
  RCP<MV> nullspace_;
};

std::string random_string( size_t length )
{
  auto randchar = []() -> char
    {
        const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}
inline bool exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


result solvekosz(int argc, char* argv[], LLA* mat) {
  result trial;
  std::string nullfile("");
  double tol = 1.0e-5;
  std::string rhsfile("");
  int maxiters = -1;
  std::string jacobiOptions="";
  std::string ilutOptions="";
  std::string belosOptions="";
  CommandLineProcessor cmdp(false,true);

  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("rhsfile",&rhsfile,"Filename for right hand side");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
  cmdp.setOption("nullfile",&nullfile,"Filename for nullspace.");
  cmdp.setOption("ilut-xml", &ilutOptions, "Specify xml file with ILUT parameters");
  cmdp.setOption("jacobi-xml",&jacobiOptions, "Specify xml file with Jacobi parameters");
  cmdp.setOption("belos-xml",&belosOptions, "Specify xml file with Jacobi parameters");
  
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    trial.success=-4;
    return trial;
  }
  
  int n=mat->n;
  Graph h(n);
  for(int i=1; i <= n; ++i) {
    h.e[i].clear();
  }
  for(int i=0; i < n; ++i) {
    LL_node *temp_node = mat->neighbor_list[i];
    for(int j=0; j < mat->num_entries[i]; ++j) {
      h.e[i+1].push_back(make_pair(temp_node->idx + 1, 1./(temp_node->weight)));
      h.e[temp_node->idx+1].push_back(make_pair(i+1,1./(temp_node->weight)));
      temp_node=temp_node->next_ptr;
    }
  }
  GraphSP g=TreeFinder::findLowStretchTree(h);
  Mat L = IO::constructMatrixFromGraph(g);
  Vec b = IO::readMMVec(rhsfile);
  Vec newx;
  KOSZ kosz(g);
  AbstractSolver S = kosz;
  int flag; FLOAT relres; int iter; vector<FLOAT> resvec; int work;
  clock_t t_start = clock();
  tie(newx, flag, relres, iter, resvec, work) = S.solve(b, 1e-12, 200);
  clock_t t_end = clock();
  h.freeMemory();
  g.freeMemory();
  kosz.freeMemory();
  L.freeMemory();
  if(flag==0) {
    trial.success=0;
  }
  else {
    trial.success=-1;
    return trial;
  }
  double tcost = static_cast<double>(t_end - t_start) / static_cast<double>(CLOCKS_PER_SEC);
    
  trial.residual=relres;
  trial.time=tcost;
  trial.flops=work;
  trial.iters=iter;
  return trial;
}

result solvetril(int argc, char *argv[], LLA* mat, std::string precond) {
  result trial;
  std::ostringstream mystream;
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  //Teuchos::RCP<Teuchos::FancyOStream> myout = Teuchos::getFancyOStream (Teuchos::rcpFromRef (mystream));
  Teuchos::RCP<std::ostream> test=Teuchos::rcpFromRef (mystream);
  GlobalMPISession mpisess(&argc,&argv,NULL);

  const ST one  = SCT::one();

  int MyPID = 0;

  
  
  Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  RCP<const Comm<int> > comm = platform.getComm();
  RCP<Node>             node = platform.getNode();

  //
  // Get test parameters from command-line processor
  //  
  bool verbose = false, proc_verbose = false, debug = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numrhs = 1;      // total number of right-hand sides to solve for
  int blocksize = 1;   // blocksize used by solver
  int blocksizefactor = 1;
  int maxiters = -1;   // maximum number of iterations for solver to use

  std::string filename("");
  std::string nullfile("");
  std::string rhsfile("");
  
  double tol = 1.0e-5;     // relative residual tolerance
  double offset = 0.01;
  int randomize = 1;
  double diagonal = 1.0;
  int forests = 1;
  
  std::string offtype = "Absolute";
  double droptol = 1e-06;
  double lof = 0.0;

  //std::string precond("none");
  std::string matrix("none");
  std::string precondform("none");

  std::string mapfile("none");

  std::string mueluOptions("mueludefault.xml");
  std::string belosOptions("belosdefault.xml");
  std::string supportgraphOptions("supportgraphdefault.xml");
  std::string ilutOptions("ilutdefault.xml");
  std::string jacobiOptions("jacobidefault.xml");
  std::string sgsOptions("sgsdefault.xml");
  std::string rilukOptions("rilukdefault.xml");
  std::string amesos2Options("amesos2default.xml");
  std::string schwarzOptions("schwarzdefault.xml");

  std::string input1("");
  std::string input2("");
  int mutatect=0;
  
  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
  cmdp.setOption("filename",&filename,"Filename for test matrix.");
  cmdp.setOption("nullfile",&nullfile,"Filename for nullspace.");
  cmdp.setOption("rhsfile",&rhsfile,"Filename for right hand side");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
  cmdp.setOption("block-size",&blocksize,"Block size to be used by the CG solver.");
  //cmdp.setOption("precond", &precond, "Preconditioner to use");
  cmdp.setOption("matrix", &matrix, "Type of matrix to use");
  cmdp.setOption("precondform", &precondform, "If the form of the preconditioner (Laplacian, etc.) is different from the matrix");
  cmdp.setOption("mapfile", &mapfile, "Specify a rowmap input file");
  cmdp.setOption("muelu-xml", &mueluOptions, "Specify xml file with muelu parameters");
  cmdp.setOption("belos-xml", &belosOptions, "Specify xml file with belos parameters");
  cmdp.setOption("supportgraph-xml", &supportgraphOptions, "Specify xml file with supportgraph parameters");
  cmdp.setOption("ilut-xml", &ilutOptions, "Specify xml file with ILUT parameters");
  cmdp.setOption("jacobi-xml",&jacobiOptions, "Specify xml file with Jacobi parameters");
  cmdp.setOption("sgs-xml",&sgsOptions, "Specify xml file with SGS parameters");
  cmdp.setOption("riluk-xml", &rilukOptions, "Specify xml file with RILUK parameters");
  cmdp.setOption("amesos2-xml", &amesos2Options, "Specify xml file with Amesos2 parameters");
  cmdp.setOption("schwarz-xml", &schwarzOptions, "Specify xml file with Additive Schwarz parameters");

  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    trial.success=-4;
    return trial;
  }
  if (debug) {
    verbose = true;
  }
  if (!verbose) {
    frequency = -1;  // reset frequency if test is not verbose
  }

  MyPID = Teuchos::rank(*comm);
  proc_verbose = ( verbose && (MyPID==0) );

  if (proc_verbose) {
    std::cout << Belos::Belos_Version() << std::endl << std::endl;
  }

  size_t num_verts = mat->n;
  std::vector<std::vector<int> > Indices (num_verts);
  std::vector<std::vector<MT> > Values (num_verts);
  Teuchos::ArrayRCP<size_t> entryct(num_verts,1);
  
  for(int i=0; i < mat->n; ++i) {
    LL_node *temp_node = mat->neighbor_list[i];
    for(int j=0; j < mat->num_entries[i]; ++j) {
      entryct[i]++;
      entryct[temp_node->idx]++;
      temp_node=temp_node->next_ptr;
    }
  }
  for(size_t row=0; row < num_verts; ++row) {
    Indices[row].resize(entryct[row]);
    Values[row].resize(entryct[row]);
  }
  for(size_t row=0; row < num_verts; ++row){
    Indices[row][0]=row;
    Values[row][0]=0;
    entryct[row]=1;
  }
  int precdigits=6;
  for(int i=0; i < mat->n; ++i) {
    LL_node *temp_node = mat->neighbor_list[i];
    for(int j=0; j < mat->num_entries[i]; ++j) {
      
      double roundv=round(pow(10.,precdigits)*temp_node->weight/pow(10.,precdigits));
      int col=temp_node->idx;
      Indices[i][entryct[i]]=col;
      Values[i][entryct[i]]=-roundv;
      Values[i][0]+=roundv;
      entryct[i]++;
      
      Indices[col][entryct[col]]=i;
      Values[col][entryct[col]]=-roundv;
      Values[col][0]+=roundv;
      entryct[col]++;
      temp_node=temp_node->next_ptr;
    }
  }

  for(size_t row=0; row < num_verts; ++row) {
    if(Values[row][0]<=0)
      {
	trial.success=-3;
	return trial;
      }
  }
  
  RCP<const Map> map = rcp(new Map(num_verts,0,comm));
  RCP<TCrsMatrix> A=rcp(new TCrsMatrix (map,map,entryct,Tpetra::StaticProfile));
  for(size_t row=0; row < num_verts; ++row) {
    Teuchos::ArrayView<int> IndicesInsert (Indices[Teuchos::as<int> (row)]);
    Teuchos::ArrayView<ST> ValuesInsert (Values[Teuchos::as<int> (row)]);
    A->insertLocalValues(row,IndicesInsert,ValuesInsert);
  }
  A->fillComplete();
  std::cout << "test4" << std::endl;
  //std::cout << "connected" << std::endl;
  //A->describe(*out,Teuchos::VERB_EXTREME);
  //print_matrix(mat);

  // Create initial vectors
  RCP<MultiVector<ST,int> > tempNullspace;
  
  if(nullfile == "ones"){
    tempNullspace=rcp(new MultiVector<ST,int>(map,1));
    MVT::MvInit(*tempNullspace,1.0);
  }
  else if(nullfile != ""){
    tempNullspace=Tpetra::MatrixMarket::Reader<MultiVector<ST,int> >::readDenseFile(nullfile,comm,map);
  }

  const RCP<MultiVector<ST,int> > Nullspace(tempNullspace);


  RCP<MultiVector<ST,int> > B, X, Xactual;
  
  if(rhsfile=="") {
    Xactual = rcp(new MultiVector<ST,int>(map,numrhs));
    X = rcp(new MultiVector<ST,int>(map,numrhs));
    MVT::MvRandom(*Xactual);
    B = rcp(new MultiVector<ST,int>(map,numrhs));
    OPT::Apply( *A, *Xactual, *B );
    //MVT::MvRandom( *X);
    MVT::MvInit(*X,0.0);
  }
  else {
    X = rcp(new MultiVector<ST,int>(map,numrhs));
    //MVT::MvRandom( *X);
    B=Tpetra::MatrixMarket::Reader<MultiVector<ST,int> >::readDenseFile(rhsfile,comm,map);
    MVT::MvInit(*X,0.0);
  }

  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->getGlobalLength();
  
  ParameterList belosList = *Teuchos::getParametersFromXmlFile(belosOptions);
  int listiters=belosList.get("Maximum Iterations", maxiters);
  

  blocksizefactor=belosList.get("Block Size Factor", blocksizefactor);  
  if(blocksizefactor!=1) {
    blocksize=NumGlobalElements/blocksizefactor;
  }
  
  belosList.set("Block Size", blocksize);

  if(maxiters!=-1 && listiters==-1){
    belosList.set("Maximum Iterations", maxiters);
  }
  else if(listiters!=-1){
    belosList.set("Maximum Iterations", listiters);
  }
  else {
    belosList.set("Maximum Iterations", NumGlobalElements/blocksize - 1);
  }

  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  int verbLevel = Belos::Errors + Belos::Warnings;
  //if (debug) {
  //  verbLevel += Belos::Debug;
  //}
  if (verbose) {
    verbLevel += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
    
  }
  verbLevel += Belos::TimingDetails + Belos::StatusTestDetails + Belos::FinalSummary;
  
  belosList.set( "Output Stream", test);
  belosList.set( "Verbosity", verbLevel );
  belosList.set( "Output Frequency",(int)maxiters/2); 

  //
  // Construct an unpreconditioned linear problem instance.
  //

  RCP<SingularOp> singularA= rcp(new SingularOp(A,Nullspace));
  Belos::LinearProblem<ST,MV,OP> problem( singularA, X, B );
  
  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    trial.success=-4;
    return trial;
  }

  //
  // ******************************************************************
  // ************Constrcut the preconditioner***********************
  // ******************************************************************
  // 

  
  typedef Ifpack2::Preconditioner<ST,int,int> prec_type;
  typedef Ifpack2::AdditiveSchwarz<Tpetra::RowMatrix<ST,int,int> > outer_prec_type;
  Teuchos::ParameterList outer;
  outer.set("order_method","rcm");
  outer.set("schwarz: combine mode", "Add");
  outer.set("schwarz: use reordering", true);
  outer.set("schwarz: overlap level", 0);
  
  //Teuchos::ParameterList inner;
  Teuchos::RCP<ParameterList > inner;
  //ParameterList outer("Ifpack2");
  //typedef Ifpack2::AdditiveSchwarz<TCrsMatrix> TPrecond;
  
  //Teuchos::RCP<inner_prec_type> prec;
  Teuchos::RCP<prec_type> prec;
  Teuchos::RCP<outer_prec_type> ASprec(new outer_prec_type(A));
  ASprec->setParameters(outer);
  if (precond != "none" && precond != "MueLu") {
    RCP<const TCrsMatrix> tempA = A;

    Ifpack2::Factory factory;
        
    std::string PrecType = "";
    //prec = Teuchos::rcp (new TPrecond (A));

    
    if(precond == "SUPPORTGRAPH" || precond == "ilut" || precond == "RILUK" || precond == "AMESOS2") {
      PrecType=precond;
      if(precond == "SUPPORTGRAPH") {
        inner = getParametersFromXmlFile(supportgraphOptions);
      }
      else if(precond == "ilut") {
        inner = Teuchos::getParametersFromXmlFile(ilutOptions);
      }
      else if(precond == "RILUK") {
        inner = Teuchos::getParametersFromXmlFile(rilukOptions);
      }
      else if(precond == "AMESOS2") {
        inner = Teuchos::getParametersFromXmlFile(amesos2Options);
        //inner->setName("Amesos2");
      }
    }
    else if(precond == "sgs" || precond == "jacobi") {
      PrecType="RELAXATION";
      if(precond == "jacobi") {
        inner = Teuchos::getParametersFromXmlFile(jacobiOptions);
      }
      if(precond == "sgs") {
        inner = Teuchos::getParametersFromXmlFile(sgsOptions);
      }
    }
    if(precond== "schwarz") {
        inner = Teuchos::getParametersFromXmlFile(schwarzOptions);
        PrecType=precond;
    }
 
    Teuchos::Time timer ("PreconditionerSetup");
    {
      Teuchos::TimeMonitor timeMon (timer, true);
      prec = factory.create(PrecType, tempA);
      prec->setParameters(*inner);
      prec->initialize();
      prec->compute();
      ASprec->setInnerPreconditioner(prec);
      ASprec->initialize();
      ASprec->compute();
    }
    
    
    problem.setLeftPrec(ASprec);
    
  }


  if(precond == "MueLu") {
    Teuchos::RCP<MueLu::TpetraOperator<ST,int,int> > mueluprec;
    Teuchos::Time timer ("PreconditionerSetup");
    {                                                         
      Teuchos::TimeMonitor timeMon (timer, true);
      //const RCP<MultiVector<ST,int> > nullspace (new MV (map,1));
      const RCP<MultiVector<ST,int> > incoords= Teuchos::null;
      //nullspace=rcp(new MultiVector<ST,int> (map, 1));
      //MVT::MvInit( *nullspace, 1.0 );
      //mueluprec = MueLu::CreateTpetraPreconditioner(A,mueluOptions,incoords,Nullspace);
    }
    std::cout << "Preconditioner computed in " << timer.totalElapsedTime();
    problem.setLeftPrec(mueluprec);
    //problem2.setLeftPrec(mueluprec);
    
  }
  

  //
  // *******************************************************************
  // *************Start the block CG iteration***********************
  // *******************************************************************
  //
  Belos::BlockCGSolMgr<ST,MV,OP> solver( rcp(&problem,false), rcp(&belosList,false) );
  //Belos::BlockCGSolMgr<ST,MV,OP> solver( rcp(&problem,false), belosList );
  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Block size used by solver: " << blocksize << std::endl;
    std::cout << "Max number of CG iterations: " << maxiters << std::endl; 
    std::cout << "Relative residual tolerance: " << tol << std::endl;
    std::cout << std::endl;
  }
  //
  // Perform solve
  //
  Belos::ReturnType ret = solver.solve();
  
  
  //
  // Compute actual residuals and error.
  //
  bool badRes = false;
  std::vector<MT> actual_resids( numrhs );
  std::vector<MT> rhs_norm( numrhs );
  MultiVector<ST,int> resid(map, numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -one, resid, one, *B, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  if (proc_verbose) {
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  }
  MT actResavg =0;
  for ( int i=0; i<numrhs; i++) {
    MT actRes = actual_resids[i]/rhs_norm[i];
    actResavg += actRes;
    if (proc_verbose) {
      std::cout<<"Problem "<<i<<" : \t"<< actRes << " " << actual_resids[i] << " " << rhs_norm[i] <<std::endl;
    }
    if (actRes > tol) badRes = true;
  }
  std::string output=mystream.str();
  //std::cout << output << std::endl << std::endl;
  int start = output.find("Operation Op*x");
  if(start==-1) {
    trial.success=-2;
    return trial;
  }
  int last = output.find(")",start);
  std::string newstring = output.substr(start,last-start);
  start=newstring.find("                    ");
  last=newstring.find("(");
  std::string time = newstring.substr(start+20,last-start-21);
  std::string iters=newstring.substr(last+1);
  //std::cout << stoi(iters)-1 << " " << time << " " << actual_resids[0]/rhs_norm[0] << std::endl;
  //std::cout << solver.getNumIters() << std::endl;
  trial.iters=solver.getNumIters();
  trial.time=stod(time);
  trial.residual=actual_resids[0]/rhs_norm[0];
  if(precond=="jacobi") {
    //Teuchos::RCP<Ifpack2::Relaxation<TCrsMatrix> > precinner = prec;
    //Ifpack2::Relaxation<TCrsMatrix> *precinner = &(*prec);
    typedef Ifpack2::Relaxation<Tpetra::RowMatrix<ST,int,int> > ourprec;
    Teuchos::RCP<ourprec> precinner = Teuchos::rcp_dynamic_cast<ourprec> (prec);
    trial.flops=A->getGlobalNumEntries()*trial.iters + precinner->getApplyFlops();
  }
  if(precond=="ilut") {
    typedef Ifpack2::ILUT<Tpetra::RowMatrix<ST,int,int> > ourprec;
    Teuchos::RCP<ourprec> precinner = Teuchos::rcp_dynamic_cast<ourprec> (prec);
    trial.flops=(precinner->getGlobalNumEntries() + A->getGlobalNumEntries()) * trial.iters;
  }
  if(precond=="sgs") {
    typedef Ifpack2::Relaxation<Tpetra::RowMatrix<ST,int,int> > ourprec;
    Teuchos::RCP<ourprec> precinner = Teuchos::rcp_dynamic_cast<ourprec> (prec);
    trial.flops=A->getGlobalNumEntries()*trial.iters + precinner->getApplyFlops();
  }
  trial.success=0;
  return trial;
}

result solve(int argc, char *argv[], LLA* mat, std::string method) {
  if(method=="kosz") {
    return solvekosz(argc,argv,mat);
  }
  else {
    return solvetril(argc,argv,mat,method);
  }
}

int main(int argc, char *argv[]) {
  std::priority_queue<result,std::vector<result>, std::greater<result> > q;
  namespace stdfs = std::experimental::filesystem;
  stdfs::path path=argv[1];
  std::string method1=argv[2];
  std::ofstream resultinfo;
  resultinfo.open(argv[3],std::ios::app);
  const stdfs::directory_iterator end{};
  std::vector<std::string> filenames;
  std::cerr << argv[1] << " " << method1 << " " << argv[3] << std::endl;
  for(stdfs::directory_iterator iter{path} ; iter != end; ++iter) {
    filenames.push_back(iter->path().string());
  }
  for(int i=0; i < filenames.size(); ++i) {
    
    LLA* matrix = MMA_read(filenames[i]);
    try{
      result trial1=solve(argc-3,argv+3,matrix,method1);
      std::cout << "finished" << std::endl;      	    
      if(trial1.success == 0) {
	
	trial1.filename=filenames[i];

	if(q.size() < 100) {
	  q.push(trial1);
	}
	else if(trial1 > q.top()){
	  result temp=q.top();
	  q.pop();
	  q.push(trial1);
	}
      }
    }
    catch(...){}
    matrix->freeMemory();
    delete matrix;
    
  }
  result temp;
  while(!q.empty()) {
    temp=q.top();
    resultinfo << temp.filename << " " << temp.flops  << std::endl;

    q.pop();
  }
    
  
  resultinfo.close();
  return 0;
}
 

