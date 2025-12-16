#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
using std::cerr;

#include "ParmParse.H"
#include "LoadBalance.H"
#include "BRMeshRefine.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "EBCellFactory.H"
#include "PoissonUtilities.H"

#include "EBFABView.H"
#include "EBDebugDump.H"

#include "EBAMRPoissonOp.H"
#include "EBLevelDataOps.H"
#include "BaseBCValue.H"
#include "BaseDomainBC.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "BaseEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "NeumannPoissonEBBC.H"

#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBLevelDataOps.H"
#include "EBSimpleSolver.H"
#include "BiCGStabSolver.H"
#include "EBEllipticLoadBalance.H"
#include "EBLoadBalance.H"
#include "LoadBalance.H"
#include "EBLevelGrid.H"
#include "CH_Timer.H"
#include "memusage.H"
#include "memtrack.H"

/********/
void solve(const PoissonParameters&  a_params)
{
  CH_TIMERS("solve driver");
  CH_TIMER("get layouts",   t1);
  CH_TIMER("define data",   t2);
  CH_TIMER("define solver", t3);
  CH_TIMER("amrmultigrid init",   t4);
  CH_TIMER("amrmultigrid solvenoinit",   t5);
  int nvar = 1;
  Vector<DisjointBoxLayout> grids;
  Vector<EBISLayout>        ebisl;

  CH_START(t1);
  getAllIrregRefinedLayouts(grids, ebisl, a_params);
  CH_STOP(t1);

  //define  data
  CH_START(t2);
  Vector<LevelData<EBCellFAB>* > phi(a_params.numLevels);
  Vector<LevelData<EBCellFAB>* > rhs(a_params.numLevels);
  Vector<LevelData<EBCellFAB>* > resid(a_params.numLevels);
  ParmParse pp;
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      EBCellFactory factory(ebisl[ilev]);
      phi[ilev] = new LevelData<EBCellFAB>(grids[ilev],nvar, a_params.ghostPhi, factory);
      rhs[ilev] = new LevelData<EBCellFAB>(grids[ilev],nvar, a_params.ghostRHS, factory);
      resid[ilev] = new LevelData<EBCellFAB>(grids[ilev],nvar, a_params.ghostRHS, factory);      

      //for now just set phi to 1 and rhs to -1.
      EBLevelDataOps::setVal(*phi[ilev], 1.0);
      EBLevelDataOps::setVal(*rhs[ilev], -1.0);
      EBLevelDataOps::setVal(*resid[ilev], 0.0);      
    }
  CH_STOP(t2);

  int relaxtype;
  pp.get("mg_relax_type", relaxtype);
  if (relaxtype !=0)
    {
      bool doLazyRelax;
      pp.get("do_lazy_relax", doLazyRelax);
      if (doLazyRelax)
        {
          pout() << "doing Lazy relaxation--one exchange for all colors" << endl;
        }
      else
        {
          pout() << "doing standard relaxation--one exchange per color" << endl;
        }
      EBAMRPoissonOp::doLazyRelax(doLazyRelax);
    }
  //create the solver
  CH_START(t3);
  AMRMultiGrid<LevelData<EBCellFAB> > solver;
  pout() << "defining  solver" << endl;

  BiCGStabSolver<LevelData<EBCellFAB> > newBottomSolver;
  newBottomSolver.m_verbosity = 0;
  Real alpha, beta;
  pp.get("alpha", alpha);
  pp.get("beta",  beta);
  defineSolver(solver, grids, ebisl, newBottomSolver, a_params, 1.e99, alpha, beta);

  CH_STOP(t3);
  pout() << "solving " << endl;
  //solve the equation
  CH_START(t4);
  solver.init(phi, rhs, a_params.maxLevel, 0);
  CH_STOP(t4);
  CH_START(t5);
  solver.solveNoInit(phi, rhs, a_params.maxLevel, 0, false);
  CH_STOP(t5);

#ifdef CH_USE_HDF5
  bool fileOut;
  pp.get("do_file_output", fileOut);
  if (fileOut)
    {
      pout() << "outputting solution to file" << endl;
      string plot_prefix("phi");
      pp.query("plot_prefix", plot_prefix);
      char charstr[100];
      sprintf(charstr, "%s.%dd.hdf5", plot_prefix.c_str(), SpaceDim);
      string filename(charstr);
      // write out phi and rhs as well as residual -- this is wrong for nvar != 1
      CH_assert(nvar == 1);
      Vector<string> names(3);
      names[0] = string("phi");
      names[1] = string("rhs");
      names[2] = string("residual");

      // compute residual
      bool homogeneousBC = false;
      // set this to be true as a reality check (can turn off later)
      bool computeNorm = true;
      int lmax = phi.size() -1;
      int lbase = 0;
      solver.computeAMRResidual(resid, phi, rhs, lmax, lbase, homogeneousBC, computeNorm);
      
      Vector<LevelData<EBCellFAB>* > plotData(phi.size());
      for (int lev=0; lev< plotData.size(); lev++)
	{
	  // 3 comps per variable: phi, rhs, residual
	  int nPlotComp = 3*nvar;
	  EBCellFactory factory(ebisl[lev]);
	  plotData[lev] = new LevelData<EBCellFAB>(grids[lev],nPlotComp, a_params.ghostPhi, factory);
	  Interval srcComps(0,nvar-1);
	  Interval phiDestComps(0,nvar-1);
	  Interval rhsDestComps(nvar,2*nvar-1);
	  Interval residDestComps(2*nvar,3*nvar-1);
	  phi[lev]->copyTo(srcComps,*plotData[lev],phiDestComps);
	  rhs[lev]->copyTo(srcComps,*plotData[lev],rhsDestComps);
	  resid[lev]->copyTo(srcComps,*plotData[lev],residDestComps);	  
	}
      
      bool replaceCovered;
      Vector<Real> coveredValues(3, 0.0);
      Real dumReal =  1.0;
      writeEBHDF5(filename, grids, plotData, names,
                  a_params.coarsestDomain.domainBox(), a_params.coarsestDx[0], dumReal, dumReal,
                  a_params.refRatio, a_params.numLevels,
                  replaceCovered, coveredValues);
      // clean up memory
      for (int lev = 0; lev < plotData.size(); lev++)
	{
	  delete plotData[lev] ;
	}
    }
#endif
  //clean up memory
  for (int ilev = 0; ilev < a_params.numLevels; ilev++)
    {
      delete phi[ilev] ;
      delete rhs[ilev] ;
      delete resid[ilev];
    }
}
/******/
int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
    CH_TIMERS("uber timers");
    CH_TIMER("define geom", t1);
    CH_TIMER("solve", t2);
#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    std::string poutBaseName = "pout";
    pp.query("poutBaseName",poutBaseName);
    setPoutBaseName(poutBaseName);
    
    PoissonParameters params;

    //read params from file
    getPoissonParameters(params);

    //define geometry from given params
    CH_START(t1);
    definePoissonGeometry(params);
    CH_STOP(t1);

    //solve the stinking problem and output everything
    CH_START(t2);
    solve(params);
    CH_STOP(t2);

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  }
  // End scoping trick

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
