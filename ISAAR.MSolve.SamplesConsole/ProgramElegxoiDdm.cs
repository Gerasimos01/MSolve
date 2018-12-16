using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers.Interfaces;

using ISAAR.MSolve.FEM; 
using ISAAR.MSolve.FEM.Elements; 
using ISAAR.MSolve.FEM.Entities; 
using ISAAR.MSolve.FEM.Materials; 
using ISAAR.MSolve.Materials; 
using ISAAR.MSolve.SamplesConsole; 
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Providers;

namespace ISAAR.MSolve.SamplesConsole
{
    public class ProgramElegxoiDdm
    {
        public static void SolveRVEExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // EPILOGH MONTELOU
            int model__builder_choice;
            model__builder_choice =31;   // 9 einai to megalo me to renumbering pou tsekaretai


            if (model__builder_choice == 1) // 
            { } //DddmExamplesBuilder.Reference1RVEExample10000(model); } //A.1 comment out osa den tous kanoume proswrina prosarmogh compatiblity
            if (model__builder_choice == 2) // 
            { }// DddmExamplesBuilder.Reference1RVEExample10000_Hexaonly(model); } //A.2
            if (model__builder_choice == 8) // 
            { RVEExamplesBuilder.Reference2RVEExample10000withRenumbering(model); }
            //if (model__builder_choice == 4) // 
            //{ DddmExamplesBuilder.Reference1RVEExample10000_Hexaonly(model); }
            if (model__builder_choice == 5) // 
            { RVEExamplesBuilder.Reference2RVEExample100_000withRenumbering_mono_hexa(model); }
            if (model__builder_choice == 6) // 
            { RVEExamplesBuilder.Reference2RVEExample500_000withRenumbering_mono_hexa(model); }
            if (model__builder_choice == 15) // 
            { RVEExamplesBuilder.Reference2RVEExample10000withRenumberingwithInput_develop2(model); }

            if (model__builder_choice == 21) // 
            { RVEExamplesBuilder.FewElementsRVECheckExample2GrapheneSheets(model); }
            if (model__builder_choice == 22) // 
            { RVEExamplesBuilder.FewElementsRVECheckExample(model); }

            if (model__builder_choice == 23) // einai to benchmark 4 
            { ParadeigmataElegxwnBuilder.HexaCantileverBuilder(model, 850); }
            if (model__builder_choice == 24) // 
            { ParadeigmataElegxwnBuilder.HexaCantileverBuilder_copyMS_222(model, 0.00219881744271988174427); }
            if (model__builder_choice == 25) // 
            { RVEExamplesBuilder.Reference2RVEExample10_000withRenumbering_mono_hexa(model); }
            if (model__builder_choice == 26) // einai to benchmark 4 
            { ParadeigmataElegxwnBuilder.HexaCantileverBuilder_material_fromMS(model, 0.00219881744271988174427); }
            if (model__builder_choice == 27) // 
            { ParadeigmataElegxwnBuilder.HexaCantileverBuilder_copyMS_222_Material_v2_DefGrad(model, 0.00219881744271988174427); }
            if (model__builder_choice == 28) // 
            { ParadeigmataElegxwnBuilder.HexaCantileverBuilder_copyMS_222_grapheneReinforced_AND_printStrainHistory(model, 0.00219881744271988174427); }


            if (model__builder_choice == 30) // einai to benchmark 4 me to neo hexa_mat
            { ParadeigmataElegxwnBuilder.Hexa_1mat_CantileverBuilder(model, 850); }
            if (model__builder_choice == 31) // Hexa8 kanoniko me NL analyzer paradeigma me Vasili Von mises
            { ParadeigmataElegxwnBuilder.HexaElementsOnlyVonMises(model); }
            if (model__builder_choice == 32) // Hexa8 kanoniko me NL analyzer paradeigma me Vasili Von mises
            { }//            IntegrationTestModelBuilders.Reference2RVEExample1000ddm_test_for_Msolve_release_version(model); } //A.3
            if (model__builder_choice == 33) // Hexa8 kanoniko me NL analyzer paradeigma me Vasili Von mises
            { ParadeigmataElegxwnBuilder.HexaElementsOnlyVonMisesHexa8NL(model); }

            bool use_domain_decomposer = false;
            if (use_domain_decomposer)
            {
                //i)
                //DddmExamplesBuilder.MakeModelDictionariesZeroBasedForDecomposer(model); //A.4

                model.ConnectDataStructures();
                // ii)
                AutomaticDomainDecomposer domainDecomposer = new AutomaticDomainDecomposer(model, 2); //2o orisma arithmoos subdomains
                domainDecomposer.UpdateModel();
            }
            else
            {
                model.ConnectDataStructures();
            }

            


            //comment section 1 palaia version
            //SolverSkyline solver = new SolverSkyline(model);
            //ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            
            ProblemStructural provider = new ProblemStructural(model, linearSystems);


            // PARADEIGMA A: LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
            //SolverSkyline2 solver = new SolverSkyline2(linearSystems[1]); //H MARIA XRHSIMOPOIEI TON sklinesolver 
            //LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            //---------------------------------------------------------------------------------------------------------------------------------

            // PARADEIGMA B: Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, 17, model.TotalDOFs);//1. increments einai to 17 (arxika eixame thesei2 26 incr)
            //PALIA DIATUPWSH: NewtonRaphsonNonLinearAnalyzer analyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystems, provider, 10, 48); 
            // NEA DIATUPWSH:
            var solver = new SolverSkyline(linearSystems[1]);
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            var increments = 2;
            var childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);
            //h epomenhgrammh einai gia paradeigma ws pros to access
            //IAnalyzer childAnalyzer2 = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);
            

            childAnalyzer.SetMaxIterations = 100;
            childAnalyzer.SetIterationsForMatrixRebuild = 1;
            //---------------------------------------------------------------------------------------------------------------------------------


            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            //childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 47 });





            //comment section 2 palaia version
            //int increments = 1;
            //Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, increments, model.TotalDOFs);//1. increments einai to 1 (arxika eixame thesei2 26 incr)
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
            //analyzer.SetMaxIterations = 100;
            //analyzer.SetIterationsForMatrixRebuild = 1;



            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


        }

        #region apla exoun ginei copy apo feat/prosthiki allagwn theoun compatibility changes A.5
        //public static void SolveDisplLoadsExample()
        //{
        //    VectorExtensions.AssignTotalAffinityCount();
        //    Model model = new Model();
        //    model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

        //    // EPILOGH MONTELOU
        //    int model__builder_choice;
        //    model__builder_choice = 33;   
            

        //    if (model__builder_choice == 33) // 
        //    { ParadeigmataElegxwnBuilder.HexaCantileverBuilderDispControl(model, 850); }
           

        //    bool use_domain_decomposer = false;
        //    if (use_domain_decomposer)
        //    {
        //        //i)
        //        DddmExamplesBuilder.MakeModelDictionariesZeroBasedForDecomposer(model);

        //        model.ConnectDataStructures();
        //        // ii)
        //        AutomaticDomainDecomposer domainDecomposer = new AutomaticDomainDecomposer(model, 2); //2o orisma arithmoos subdomains
        //        domainDecomposer.UpdateModel();
        //    }
        //    else
        //    {
        //        model.ConnectDataStructures();
        //    }

        //    Dictionary<int, Vector> uInitialFreeDOFDisplacementsPerSubdomain = new Dictionary<int, Vector>();
        //    uInitialFreeDOFDisplacementsPerSubdomain.Add(model.SubdomainsDictionary[1].ID, new Vector(model.SubdomainsDictionary[1].TotalDOFs));// prosoxh sto Id twn subdomain
        //    Dictionary<int, Node> boundaryNodes = new Dictionary<int, Node>();
        //    for (int k = 17; k < 21; k++)
        //    {
        //        boundaryNodes.Add(model.NodesDictionary[k].ID, model.NodesDictionary[k]);
        //    }
        //    Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements=new Dictionary<int, Dictionary<DOFType, double>>();
        //    Dictionary<DOFType, double> initialConvergedBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
        //    initialConvergedBoundaryNodalDisplacements.Add(DOFType.X, 0);
        //    for (int k = 17; k < 21; k++)
        //    {
        //        initialConvergedBoundaryDisplacements.Add(model.NodesDictionary[k].ID, initialConvergedBoundaryNodalDisplacements);
        //    }
        //    Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
        //    double[] prescribedDisplacmentXValues = new double[4] { 7.81614E-01, 7.07355E-01, 7.81614E-01, 7.07355E-01 };
        //    for (int k = 17; k < 21; k++)
        //    {
        //        Dictionary<DOFType, double> totalBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
        //        totalBoundaryNodalDisplacements.Add(DOFType.X, prescribedDisplacmentXValues[k-17]);
        //        totalBoundaryDisplacements.Add(model.NodesDictionary[k].ID, totalBoundaryNodalDisplacements);
        //    }

        //    ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();
        //    Dictionary<int, EquivalentContributionsAssebler> equivalentContributionsAssemblers = new Dictionary<int, EquivalentContributionsAssebler>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
        //    equivalentContributionsAssemblers.Add(model.SubdomainsDictionary[1].ID, new  EquivalentContributionsAssebler(model.SubdomainsDictionary[1], elementProvider));

        //    //comment section 1 palaia version
        //    //SolverSkyline solver = new SolverSkyline(model);
        //    //ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);

        //    var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
        //    linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces); // elegxos me model.subdomainsDictionary[1]

        //    ProblemStructural provider = new ProblemStructural(model, linearSystems);


        //    // PARADEIGMA A: LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
        //    //SolverSkyline2 solver = new SolverSkyline2(linearSystems[1]); //H MARIA XRHSIMOPOIEI TON sklinesolver 
        //    //LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
        //    //---------------------------------------------------------------------------------------------------------------------------------

        //    // PARADEIGMA B: Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, 17, model.TotalDOFs);//1. increments einai to 17 (arxika eixame thesei2 26 incr)
        //    //PALIA DIATUPWSH: NewtonRaphsonNonLinearAnalyzer analyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystems, provider, 10, 48); 
        //    // NEA DIATUPWSH:
        //    var solver = new SolverSkyline(linearSystems[1]);
        //    var linearSystemsArray = new[] { linearSystems[1] };
        //    var subdomainUpdaters = new[] { new NonLinearSubdomainUpdaterWithInitialConditions(model.Subdomains[0]) };
        //    var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

        //    var increments = 2;
        //    var childAnalyzer = new NewtonRaphsonNonLinearAnalyzerDevelopCopy(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs, uInitialFreeDOFDisplacementsPerSubdomain,
        //        boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, equivalentContributionsAssemblers);
        //    //h epomenhgrammh einai gia paradeigma ws pros to access
        //    //IAnalyzer childAnalyzer2 = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);


        //    childAnalyzer.SetMaxIterations = 100;
        //    childAnalyzer.SetIterationsForMatrixRebuild = 1;
        //    //---------------------------------------------------------------------------------------------------------------------------------


        //    StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

        //    //childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 47 });





        //    //comment section 2 palaia version
        //    //int increments = 1;
        //    //Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, increments, model.TotalDOFs);//1. increments einai to 1 (arxika eixame thesei2 26 incr)
        //    //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
        //    //analyzer.SetMaxIterations = 100;
        //    //analyzer.SetIterationsForMatrixRebuild = 1;



        //    parentAnalyzer.BuildMatrices();
        //    parentAnalyzer.Initialize();
        //    parentAnalyzer.Solve();


        //}

        //public static void SolveDisplLoadsExampleRestartAnalysis()
        //{
        //    VectorExtensions.AssignTotalAffinityCount();
        //    Model model = new Model();
        //    model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

        //    // EPILOGH MONTELOU
        //    int model__builder_choice;
        //    model__builder_choice = 33;
        //    var increments = 1;
        //    double displFactor = 0.5;
        //    double displFactor2 = 1;

        //    if (model__builder_choice == 33) // 
        //    { ParadeigmataElegxwnBuilder.HexaCantileverBuilderDispControl(model, 850); }


        //    bool use_domain_decomposer = false;
        //    if (use_domain_decomposer)
        //    {
        //        //i)
        //        DddmExamplesBuilder.MakeModelDictionariesZeroBasedForDecomposer(model);

        //        model.ConnectDataStructures();
        //        // ii)
        //        AutomaticDomainDecomposer domainDecomposer = new AutomaticDomainDecomposer(model, 2); //2o orisma arithmoos subdomains
        //        domainDecomposer.UpdateModel();
        //    }
        //    else
        //    {
        //        model.ConnectDataStructures();
        //    }

        //    #region Rve Free Dofs DIsplacement Dictionary Creation 1st Displacement Load
        //    // u (or uplusDu) initial 
        //    Dictionary<int, Vector> uInitialFreeDOFDisplacementsPerSubdomain = new Dictionary<int, Vector>();
        //    uInitialFreeDOFDisplacementsPerSubdomain.Add(model.SubdomainsDictionary[1].ID, new Vector(model.SubdomainsDictionary[1].TotalDOFs));// prosoxh sto Id twn subdomain
        //    Dictionary<int, Node> boundaryNodes = new Dictionary<int, Node>();
        //    for (int k = 17; k < 21; k++)
        //    {
        //        boundaryNodes.Add(model.NodesDictionary[k].ID, model.NodesDictionary[k]);
        //    }
        //    Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
        //    Dictionary<DOFType, double> initialConvergedBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
        //    initialConvergedBoundaryNodalDisplacements.Add(DOFType.X, 0);
        //    for (int k = 17; k < 21; k++)
        //    {
        //        initialConvergedBoundaryDisplacements.Add(model.NodesDictionary[k].ID, initialConvergedBoundaryNodalDisplacements);
        //    }
        //    Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
        //    double[] prescribedDisplacmentXValues = new double[4] { 7.81614E-01, 7.07355E-01, 7.81614E-01, 7.07355E-01 };
        //    for (int k = 17; k < 21; k++)
        //    {
        //        Dictionary<DOFType, double> totalBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
        //        totalBoundaryNodalDisplacements.Add(DOFType.X, displFactor*prescribedDisplacmentXValues[k - 17]);
        //        totalPrescribedBoundaryDisplacements.Add(model.NodesDictionary[k].ID, totalBoundaryNodalDisplacements);
        //    }
        //    #endregion

        //    #region Creation of nessesary analyzers and solution 
        //    ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();
        //    Dictionary<int, EquivalentContributionsAssebler> equivalentContributionsAssemblers = new Dictionary<int, EquivalentContributionsAssebler>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
        //    equivalentContributionsAssemblers.Add(model.SubdomainsDictionary[1].ID, new EquivalentContributionsAssebler(model.SubdomainsDictionary[1], elementProvider));
           
        //    var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
        //    linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces); // Comment MS: Gia to multiscale mporoume apla na tou perasoume ena keno dianusma me diastash force afou sto MS  den uparxei epivolh
        //    // fortiou mesw loads, to forces to ftiaxnei to model.connectDataStructures elegxos me model.subdomainsDictionary[1]

        //    ProblemStructural provider = new ProblemStructural(model, linearSystems);
        //    var solver = new SolverSkyline(linearSystems[1]);
        //    var linearSystemsArray = new[] { linearSystems[1] };
        //    var subdomainUpdaters = new[] { new NonLinearSubdomainUpdaterWithInitialConditions(model.Subdomains[0]) };
        //    var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

        //    var childAnalyzer = new NewtonRaphsonNonLinearAnalyzerDevelopCopy(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs, uInitialFreeDOFDisplacementsPerSubdomain,
        //        boundaryNodes, initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, equivalentContributionsAssemblers);            
        //    childAnalyzer.SetMaxIterations = 100;
        //    childAnalyzer.SetIterationsForMatrixRebuild = 1;
            

        //    StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);            
        //    parentAnalyzer.BuildMatrices();
        //    parentAnalyzer.Initialize();
        //    parentAnalyzer.Solve();
        //    #endregion

        //    #region Rve Free Dofs DIsplacement Dictionary Creation 2nd Displacement Load
        //    // u (or uplusDu) initial 
        //    uInitialFreeDOFDisplacementsPerSubdomain = childAnalyzer.GetConvergedSolutionVectorOfFreeDofs();

        //    initialConvergedBoundaryDisplacements = totalPrescribedBoundaryDisplacements;


        //    totalPrescribedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
        //    for (int k = 17; k < 21; k++)
        //    {
        //        Dictionary<DOFType, double> totalBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
        //        totalBoundaryNodalDisplacements.Add(DOFType.X, displFactor2 * prescribedDisplacmentXValues[k - 17]);
        //        totalPrescribedBoundaryDisplacements.Add(model.NodesDictionary[k].ID, totalBoundaryNodalDisplacements);
        //    }
        //    #endregion

        //    #region Creation of nessesary analyzers and solution 
        //    linearSystems[1].RHS.Clear();
        //    //linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces); // Comment MS: Gia to multiscale mporoume apla na tou perasoume ena keno dianusma me diastash force afou sto MS  den uparxei epivolh
        //    // fortiou mesw loads, to forces to ftiaxnei to model.connectDataStructures elegxos me model.subdomainsDictionary[1]

        //    childAnalyzer = new NewtonRaphsonNonLinearAnalyzerDevelopCopy(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs, uInitialFreeDOFDisplacementsPerSubdomain,
        //        boundaryNodes, initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, equivalentContributionsAssemblers);
        //    childAnalyzer.SetMaxIterations = 100;
        //    childAnalyzer.SetIterationsForMatrixRebuild = 1;


        //    parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
        //    parentAnalyzer.BuildMatrices();
        //    parentAnalyzer.Initialize();
        //    parentAnalyzer.Solve();
        //    #endregion



        //}
        #endregion

        //static void Main(string[] args)
        //{
        //    SolveRVEExample(); //|
        //}

    }
}
