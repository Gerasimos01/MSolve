using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.LinearSystems;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Output;
using ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses;
using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.SamplesConsole
{
    public class SeparateCodeCheckingClass7_b_b
    {
        //origin: epilush paradeigmatwn antistoixwn ths SeparateCodeCheckingClass5b_b dld me RveGrShMultipleSeparatedDevelopb
        // OneRveExample.Check_Graphene_rve_serial() tha xrhsimopoiithei gia ton orismo tou provlhmatos

        //edw xreiazetai rve apo ton builder xwris connectdata structures kai  xwris defineAppropriateConstraints
        public static void NLRVEStrainParralelSolution()
        {
            var rveBuilder = new RveGrShMultipleSeparatedDevelopb(1, true);
            double[] uc;

            var microstructure = new MicrostructureDefGrad3D(rveBuilder, false, 1);


        }

        
    }

    public class SeparateCodeCheckingClass_c_alte_develop
    {
        //origin: epilush paradeigmatwn antistoixwn ths SeparateCodeCheckingClass5b_b dld me RveGrShMultipleSeparatedDevelopb
        // OneRveExample.Check_Graphene_rve_serial() tha xrhsimopoiithei gia ton orismo tou provlhmatos

        //edw xreiazetai rve apo ton builder xwris connectdata structures kai  xwris defineAppropriateConstraints
        public static void NLRVEStrainParralelSolution()
        {
            //parallel
            var rveBuilder2 = new RveGrShMultipleSeparated_c_alteDevelop(1, true);
            
            var microstructure2 = new MicrostructureDefGrad3Ddevelop(rveBuilder2, false, 1);
            double[] DGvec2 = new double[9] { 1.00000000001, 1, 1, 0, 0, 0, 0, 0, 0 };

            microstructure2.UpdateMaterialDevelop(DGvec2);

            double[] solution2 = microstructure2.iterSolutionGathered[0].CopyToArray();
            var solutions2 = microstructure2.iterSolutions[0][0].CopyToArray();

            var cornerNodesIds = rveBuilder2.CornerNodesIds;
            double[] cornerNodeSolution2 = new double[3 * cornerNodesIds.Keys.Count];
            int thesi = 0;
            foreach(int cornerNodeID in cornerNodesIds.Keys)
            {
                var cornerNode = microstructure2.model.NodesDictionary[cornerNodeID];
                var dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
                foreach(StructuralDof dof in dofs)
                {
                    int globalDof = microstructure2.model.GlobalDofOrdering.GlobalFreeDofs[cornerNode, dof];
                    double u = solution2[globalDof];
                    cornerNodeSolution2[thesi] = u;
                    thesi++;
                }

            }

            // serial 
            var rveBuilder = new RveGrShMultipleSeparated_c_alteDevelop(1, false);

            var microstructure = new MicrostructureDefGrad3Ddevelop(rveBuilder, false, 1);
            double[] DGvec = new double[9] { 1.00000000001, 1, 1, 0, 0, 0, 0, 0, 0 };

            microstructure.UpdateMaterialDevelop(DGvec);

            var solution = microstructure.iterSolutions[0][1].CopyToArray();

            double[] cornerNodeSolution = new double[3 * cornerNodesIds.Keys.Count];
            thesi = 0;
            foreach (int cornerNodeID in cornerNodesIds.Keys)
            {
                var cornerNode = microstructure.model.NodesDictionary[cornerNodeID];
                var dofs = new StructuralDof[3] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
                foreach (StructuralDof dof in dofs)
                {
                    int globalDof = microstructure.model.GlobalDofOrdering.GlobalFreeDofs[cornerNode, dof];
                    double u = solution[globalDof];
                    cornerNodeSolution[thesi] = u;
                    thesi++;
                }

            }

            //microstructure.UpdateMaterialDevelop(new double[9] { 1, 1, 1, 0, 0, 0, 0, 0, 0 });

        }

        public static void NLRVEStrainParralelSolution2()
        {
            var rveBuilder = new RveGrShMultipleSeparated_c_alteDevelop(1, false);
            double[] uc;

            var microstructure = new MicrostructureDefGrad3Ddevelop(rveBuilder, false, 1);
            double[] DGvec = new double[9] { 1.00001, 1, 1, 0, 0, 0, 0, 0, 0 };

            microstructure.UpdateMaterialDevelop(DGvec);

            double[] solution = microstructure.iterSolutionGathered[0].CopyToArray();
            var solutions = microstructure.iterSolutions[0][0].CopyToArray();

            microstructure.UpdateMaterialDevelop(new double[9] { 1, 1, 1, 0, 0, 0, 0, 0, 0 });

            var model = microstructure.model;

            rveBuilder.ReassignModel(model, true);
            microstructure.createSolver = rveBuilder.GetParallelSolver;

            microstructure.UpdateMaterialDevelop(DGvec);

            double[] solution2 = microstructure.iterSolutionGathered[0].CopyToArray();
            var solutions2 = microstructure.iterSolutions[0][0].CopyToArray();
        }
    }
}

