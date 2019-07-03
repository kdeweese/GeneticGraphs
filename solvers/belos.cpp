//#include <quadmath.h>
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
//typedef Kokkos::DefaultNode::DefaultNodeType Node;
//typedef Kokkos::SerialNode Node;
typedef Tpetra::CrsMatrix<ST,int,int> TCrsMatrix;
//typedef Zoltan2::XpetraCrsMatrixAdapter<TCrsMatrix> SparseMatrixAdapter;

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

int main(int argc, char *argv[]) {

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

  std::string precond("none");
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
  cmdp.setOption("precond", &precond, "Preconditioner to use");
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
    return -1;
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

  RCP<const Map> map;
  RCP<const Map> colMap;
  RCP<TCrsMatrix> A;
 

  if(mapfile != "none"){
    map = Tpetra::MatrixMarket::Reader<TCrsMatrix>::readMapFile(mapfile,comm,node,false,false);
    A = Tpetra::MatrixMarket::Reader<TCrsMatrix>::readSparseFile(filename,map,colMap,map,map,true,false,false);
  }

  else{
    A = Tpetra::MatrixMarket::Reader<TCrsMatrix>::readSparseFile(filename,comm,node,true,false,false);
    map = A->getRowMap();
  }

  
  
  RCP<MultiVector<ST,int> > tempNullspace;
  
  if(nullfile == "ones"){
    tempNullspace=rcp(new MultiVector<ST,int>(map,1));
    MVT::MvInit(*tempNullspace,1.0);
  }
  else if(nullfile != ""){
    tempNullspace=Tpetra::MatrixMarket::Reader<MultiVector<ST,int> >::readDenseFile(nullfile,comm,map);
  }

  const RCP<MultiVector<ST,int> > Nullspace(tempNullspace);

  size_t num_verts=A->getNodeNumRows();
  size_t num_entries;
  size_t max_num_entries = A->getNodeMaxNumRowEntries();

  std::vector<ST> valuestemp (max_num_entries);
  std::vector<int> indicestemp (max_num_entries);

  Teuchos::ArrayView<ST> values (valuestemp);
  Teuchos::ArrayView<int> indices (indicestemp);

  std::vector<std::vector<int> > Indices (num_verts);
  std::vector<std::vector<MT> > Values (num_verts);

  Teuchos::ArrayRCP<size_t> entryct(num_verts,1);
  
  for(size_t row=0; row < num_verts; ++row) {
    A->getLocalRowCopy(row,indices,values,num_entries);
    for(size_t colIndex=0; colIndex < num_entries;++colIndex) {
      entryct[row]++;
    }
  }
  for(size_t row=0; row < num_verts; ++row) {
    Indices[row].resize(entryct[row]);
    Values[row].resize(entryct[row]);
  }
  for(size_t row=0; row < num_verts; ++row){
    Indices[row][0]=row;
    Values[row][0]=0;
  }
  int precdigits=6;
  for(size_t row=0; row< num_verts; ++row) {
    A->getLocalRowCopy(row,indices,values,num_entries);
    for(size_t colIndex=0; colIndex < num_entries; ++colIndex) {
      double roundv=round(pow(10.,precdigits)*values[colIndex]/pow(10.,precdigits));
      Indices[row][colIndex+1]=indices[colIndex];
      Values[row][colIndex+1]=-roundv;
      Values[row][0]+=roundv;
    }
  }
  for(size_t row=0; row < num_verts; ++row) {
    if(Values[row][0]<=0)
      {
	std::cout << -3 << " " << -3 << " " << -3 << std::endl;
	return -1;
      }
  }

  RCP<TCrsMatrix> newA=rcp(new TCrsMatrix (A->getRowMap(),A->getColMap(),entryct,Tpetra::StaticProfile));
  for(size_t row=0; row < num_verts; ++row) {
    Teuchos::ArrayView<int> IndicesInsert (Indices[Teuchos::as<int> (row)]);
    Teuchos::ArrayView<ST> ValuesInsert (Values[Teuchos::as<int> (row)]);
    newA->insertLocalValues(row,IndicesInsert,ValuesInsert);
  }
  newA->fillComplete();
  A=newA;

  // Create initial vectors
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
  //if(maxiters==-1) {
    //belosList.set("Maximum Iterations", NumGlobalElements/blocksize - 1);
  //
  
  //belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
  //belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
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
  /*if (verbose) {
    if (frequency > 0) {
      belosList.set( "Output Frequency", frequency );
    }
    }*/
  //
  // Construct an unpreconditioned linear problem instance.
  //
  
  RCP<SingularOp> singularA= rcp(new SingularOp(A,Nullspace));
  Belos::LinearProblem<ST,MV,OP> problem( singularA, X, B );
  
  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
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
  outer.set("schwarz: use reordering", false);
  outer.set("schwarz: overlap level", 0);
  
  Teuchos::RCP<ParameterList > inner;
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
        //Teuchos::RCP<Ifpack2::Relaxation<TCrsMatrix> > relaxation_ptr = Teuchos::rcp( new Ifpack2::Relaxation<TCrsMatrix> (A)); 
        //prec->setInnerPreconditioner(relaxation_ptr);
      
    }
    if(precond== "SCHWARZ") {
        inner = Teuchos::getParametersFromXmlFile(schwarzOptions);
        PrecType=precond;
    }
 
    //outer.set("inner preconditioner parameters", *inner);
    //prec->setParameters(outer);
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
  
  //Pre solve (test for stagnation)
  /*ParameterList belosList2(belosList);
  int actualiters=belosList.get("Maximum Iterations", maxiters);
  belosList2.set("Maximum Iterations", actualiters-1);
  belosList2.set( "Output Stream", Teuchos::rcp(&std::cout,false));
  belosList2.set( "Verbosity", 0 );
  Belos::BlockCGSolMgr<ST,MV,OP> solver2( rcp(&problem2,false), rcp(&belosList2,false));
  Belos::ReturnType ret2 = solver2.solve();
  std::vector<MT> actual_resids2( numrhs );
  std::vector<MT> rhs_norm2( numrhs );
  MultiVector<ST,int> resid2(map, numrhs);
  OPT::Apply( *A, *preX, resid2 );
  MVT::MvAddMv( -one, resid2, one, *preB, resid2 );
  MVT::MvNorm( resid2, actual_resids2 );
  MVT::MvNorm( *preB, rhs_norm2 );
  */


  //solver.describe(*out,Teuchos::VERB_EXTREME);
  //Teuchos::Array< Teuchos::RCP < Teuchos::Time > > timers = solver.getTimers();
  //std::cout << "Elapsed Time " << (*(timers.begin()))->totalElapsedTime() << std::endl;
  //std::ofstream solutionoutput;
  //solutionoutput.open("difference.txt");
  //Teuchos::RCP<Teuchos::FancyOStream> out2 = Teuchos::getFancyOStream (Teuchos::rcpFromRef (solutionoutput));
  //Xactual->update(-1,*X,1);
  //Xactual->elementWiseMultiply(1,*(XactualRecip->getVector(0)),*Xactual,0);
  //Xactual->abs(*Xactual);
  //Xactual->describe(*out2,Teuchos::VERB_EXTREME);
  //solutionoutput.close();
  
  //
  // Compute actual residuals and error.
  //
  bool badRes = false;
  std::vector<MT> actual_resids( numrhs );
  std::vector<MT> rhs_norm( numrhs );
  MultiVector<ST,int> resid(map, numrhs);
  OPT::Apply( *singularA, *X, resid );
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
    std::cout << -2 << " " << -2 << " " << -2 << std::endl;
    return -1;
  }
  int last = output.find(")",start);
  std::string newstring = output.substr(start,last-start);
  start=newstring.find("                    ");
  last=newstring.find("(");
  std::string time = newstring.substr(start+20,last-start-21);
  std::string iters=newstring.substr(last+1);
  //if(actual_resids[0]/rhs_norm[0] < actual_resids2[0]/rhs_norm2[0] && actual_resids[0]/rhs_norm[0] < 1) 
  std::cout << stoi(iters)-1 << " " << time << " " << actual_resids[0]/rhs_norm[0] << std::endl;
  //else
  //std::cout << -1 << " " << -1 << " " << -1 << std::endl;
  //return 0;
  if (ret!=Belos::Converged || badRes) {
    if (proc_verbose) {
      std::cout << "\nEnd Result: TEST FAILED" << std::endl;	
    }
    return -1;
  }
  


  //
  // Default return value
  //
  if (proc_verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
  std::vector<MT> actual_errors( numrhs );
  std::vector<MT> lhs_norm( numrhs );
  MultiVector<ST,int> error(map, numrhs);
  MVT::MvAddMv( -one, *X, one, *Xactual, error );
  MVT::MvNorm( error, actual_errors );
  MVT::MvNorm( *Xactual, lhs_norm );
  if (proc_verbose) {
    std::cout<< "---------- Actual Errors (normalized) ----------"<<std::endl<<std::endl;
  }
  MT actErravg=0;
  for ( int i=0; i<numrhs; i++) {
    MT actErr = actual_errors[i]/lhs_norm[i];
    actErravg+=actErr;
    if (proc_verbose) {
      std::cout<<"Problem "<<i<<" : \t"<< actErr <<std::endl;
    }
  }
  std::cout << "Actual error average of system " << actErravg/numrhs << std::endl;
  
  if(nullfile != "" && nullfile!="ones") {
    RCP<MultiVector<ST,int> > NullspaceRecip = rcp(new MultiVector<ST,int>(map,1));
    NullspaceRecip->reciprocal(*Nullspace);
    X->elementWiseMultiply(1.0,*(NullspaceRecip->getVector(0)),*X,0.0);
    Xactual->elementWiseMultiply(1.0,*(NullspaceRecip->getVector(0)),*Xactual,0.0);
    MVT::MvAddMv( -one, *X, one, *Xactual, error );
    MVT::MvNorm( error, actual_errors );
    MVT::MvNorm( *Xactual, lhs_norm );
 
    if (proc_verbose) {
      std::cout<< "---------- Actual Errors of original Laplacian system (normalized) ----------"<<std::endl<<std::endl;
    }
    actErravg=0;
    for ( int i=0; i<numrhs; i++) {
      MT actErr = actual_errors[i]/lhs_norm[i];
      actErravg+=actErr;
      if (proc_verbose) {
	std::cout<<"Problem "<<i<<" : \t"<< actErr <<std::endl;
      }
    } 
    std::cout << "Actual error average of original Laplacian system " << actErravg/numrhs << std::endl;
  }


  if(precond != "none" && precond != "MueLu")
    {
      prec->describe(*out,Teuchos::VERB_HIGH);
    }
  
  return 0;
  //
} // end test_bl_cg_hb.cpp

