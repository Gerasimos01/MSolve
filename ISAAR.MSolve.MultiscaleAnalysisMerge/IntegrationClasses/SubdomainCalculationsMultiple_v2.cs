﻿using System;
using System.Collections.Generic;
//using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
//using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.Solvers.Skyline;
using ISAAR.MSolve.Solvers.Interfaces;
using System.Linq;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Commons;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public static class SubdomainCalculationsMultiple_v2
    {
        //TODO erase this one
        public static Dictionary<int, double[][]> CalculateKffinverseKfpDqSubdomains(Dictionary<int,double[][]> KfpDqSubdomains , Model model, IElementMatrixProvider elementProvider, 
            IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes, ISolver solver, Dictionary<int, ILinearSystem> linearSystems)
        {
            #region Creation of solution vectors structure
                Dictionary<int, double[][]> f2_vectorsSubdomains = new Dictionary<int, double[][]>();
            foreach(int subdomainID in KfpDqSubdomains.Keys)
            {
                f2_vectorsSubdomains.Add(subdomainID, new double[KfpDqSubdomains[subdomainID].GetLength(0)][]);
            }
            #endregion

            #region Creation of linear systems with no RHS (first RHS value can be assigned too )
            ILinearSystem[] seclinearSystems = new ILinearSystem[linearSystems.Count];
            int counter = 0;
            foreach (ILinearSystem subdomain in linearSystems.Values)
            {
                seclinearSystems[counter] = new SkylineLinearSystem(subdomain.ID, new double[KfpDqSubdomains[subdomain.ID][0].GetLength(0)]);
                seclinearSystems[counter].Matrix = subdomain.Matrix;
            }
            #endregion

            #region creation of solver
            var secSolver = new SolverSkyline(seclinearSystems[0]);
            secSolver.Initialize();
            //int subdomainID = 1;
            //for (int k=0; k< KfpDqSubdomains[subdomainID].GetLength(0); k++)
            //{
            //    secSolver[k] = new SolverSkyline(seclinearSystems[k]);
            //    secSolver[k].Initialize();
            //    secSolver[k].Solve();
            //}
            #endregion

            #region Consecutively(for macroscaleVariableDimension times) Set proper right hand side. Solve. Copy solution in output vector 
            int oneSubomainID = seclinearSystems[0].ID;
            for (int k = 0; k < scaleTransitions.MacroscaleVariableDimension(); k++) //KfpDqSubdomains[linearSystems[0].ID].GetLength(0)=Mac
            {
                #region Set proper RHS 
                //var globalRHS = new Vector(model.TotalDOFs); //TODO: uncoomment if globalRHS is needed for solver
                foreach (ILinearSystem secSubdomain in seclinearSystems)
                {
                    var subdomainRHS = ((Numerical.LinearAlgebra.Vector)secSubdomain.RHS);
                    subdomainRHS.Clear();
                    subdomainRHS.Add( new Numerical.LinearAlgebra.Vector(KfpDqSubdomains[secSubdomain.ID][k]));
                    //mappings[seclinearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == secSubdomain.ID).Index].SubdomainToGlobalVector(subdomainRHS.Data, globalRHS.Data);
                    //TODO: uncoomment if globalRHS is needed for solver
                }
                #endregion

                #region Solve
                secSolver.Solve();
                #endregion

                #region Copy solution in output vector
                foreach (ILinearSystem secSubdomain in seclinearSystems)
                {
                    f2_vectorsSubdomains[secSubdomain.ID][k] = new double[secSubdomain.Solution.Length];
                    secSubdomain.Solution.CopyTo(f2_vectorsSubdomains[secSubdomain.ID][k], 0);                    
                }
                #endregion
            }
            #endregion

            return f2_vectorsSubdomains;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="KfpDqSubdomains"></param>
        /// <param name="model"></param>
        /// <param name="elementProvider"></param>
        /// <param name="scaleTransitions"></param>
        /// <param name="boundaryNodes"></param>
        /// <param name="solver">
        /// <paramref name="solver"/>.<see cref="ISolver_v2.Initialize"/> must already have been called. Also the linear system matrices must already have been set.
        /// </param>
        

        public static Dictionary<int, double[][]>  CalculateKfreeprescribedDqMultiplicationsMultiple(Model model, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, double[][]> KfpDqSubdomains = new Dictionary<int, double[][]>();

            foreach (Subdomain subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)" tou opoiu ta stoixeia ginontai access kai me ID
            {
                KfpDqSubdomains.Add(subdomain.ID, SubdomainCalculations_v2.CalculateKfreeprescribedDqMultiplications(subdomain, elementProvider, scaleTransitions, boundaryNodes));                    
            }

            return KfpDqSubdomains;
        }

        

        public static Dictionary<int, double[][]> CalculateKppDqSubdomainsMultiplications(Model model, IElementMatrixProvider elementProvider, IScaleTransitions scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, double[][]> KppDqVectorsSubdomains = new Dictionary<int, double[][]>();

            foreach (Subdomain subdomain in model.Subdomains)
            {
                KppDqVectorsSubdomains.Add(subdomain.ID, SubdomainCalculations_v2.CalculateKppDqMultiplications(subdomain, elementProvider, scaleTransitions, boundaryNodes));
            }

            return KppDqVectorsSubdomains;
        }

        public static double[][] CombineMultipleSubdomainsIntegrationVectorsIntoTotal(Dictionary<int, double[][]> VectorsSubdomains, IScaleTransitions scaleTransitions)
        {


            double[][] totalVectors = new double[scaleTransitions.MacroscaleVariableDimension()][]; //or VectorsSubdomains.getLength(0);
            int oneSubdomainID = VectorsSubdomains.Keys.ElementAt(0);
            for (int i1 = 0; i1 < scaleTransitions.MacroscaleVariableDimension(); i1++)
            {
                totalVectors[i1] = new double[VectorsSubdomains[oneSubdomainID][i1].GetLength(0)];
            }


            foreach (int subdomainID in VectorsSubdomains.Keys)
            {
                for (int i1 = 0; i1 < scaleTransitions.MacroscaleVariableDimension(); i1++)
                {
                    for (int i2 = 0; i2 < VectorsSubdomains[subdomainID][i1].GetLength(0); i2++)
                    {
                        totalVectors[i1][i2] += VectorsSubdomains[subdomainID][i1][i2];
                    }
                }
            }

            return totalVectors;

        }       
     

        #region v2 methods
        public static Dictionary<int, double[]> CalculateFppReactionsVectorSubdomains_v2(Model_v2 model, IElementMatrixProvider elementProvider,
            IScaleTransitions_v2 scaleTransitions, Dictionary<int, Node> boundaryNodes, Dictionary<int, IVectorView> solution, Dictionary<int, IVectorView> dSolution,
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {
            Dictionary<int, double[]> FppReactionVectorSubdomains = new Dictionary<int, double[]>();

            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                FppReactionVectorSubdomains.Add(subdomain.ID, SubdomainCalculations_v2.CalculateFppReactionsVector_v2(subdomain, elementProvider, scaleTransitions, boundaryNodes,
                solution[subdomain.ID], dSolution[subdomain.ID], initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements));
            }

            return FppReactionVectorSubdomains;

        }

        internal static double[] CombineMultipleSubdomainsStressesIntegrationVectorsIntoTotal_v2(Dictionary<int, double[]> fppReactionVectorSubdomains)
        {

            double[] totalVector = new double[fppReactionVectorSubdomains.ElementAt(0).Value.GetLength(0)];

            foreach (int subdomainID in fppReactionVectorSubdomains.Keys)
            {
                for (int i1 = 0; i1 < totalVector.GetLength(0); i1++)
                {
                    totalVector[i1] += fppReactionVectorSubdomains[subdomainID][i1];
                }
            }

            return totalVector;

        }

        public static Dictionary<int, double[][]> CalculateKffinverseKfpDqSubdomains_v2(Dictionary<int, double[][]> KfpDqSubdomains, Model_v2 model, IElementMatrixProvider elementProvider,
            IScaleTransitions_v2 scaleTransitions, Dictionary<int, Node> boundaryNodes, ISolver_v2 solver)
        {
            IReadOnlyDictionary<int, ILinearSystem_v2> linearSystems = solver.LinearSystems;

            #region Creation of solution vectors structure
            Dictionary<int, double[][]> f2_vectorsSubdomains = new Dictionary<int, double[][]>();
            foreach (int subdomainID in KfpDqSubdomains.Keys)
            {
                f2_vectorsSubdomains.Add(subdomainID, new double[KfpDqSubdomains[subdomainID].GetLength(0)][]);
            }
            #endregion

            //#region Creation of linear systems with no RHS (first RHS value can be assigned too )
            //ILinearSystem_v2[] seclinearSystems = new ILinearSystem_v2[linearSystems.Count];
            //int counter = 0;
            //foreach (ILinearSystem_v2 subdomain in linearSystems.Values)
            //{
            //    seclinearSystems[counter] = new SkylineLinearSystem(subdomain.ID, new double[KfpDqSubdomains[subdomain.ID][0].GetLength(0)]);
            //    seclinearSystems[counter].Matrix = subdomain.Matrix;
            //}
            //#endregion

            //#region creation of solver
            //var secSolver = new SolverSkyline(seclinearSystems[0]);
            //secSolver.Initialize();

            //#endregion

            #region Consecutively(for macroscaleVariableDimension times) Set proper right hand side. Solve. Copy solution in output vector 
            int oneSubomainID = linearSystems.First().Value.Subdomain.ID;           //seclinearSystems[0].ID;
            for (int k = 0; k < scaleTransitions.MacroscaleVariableDimension(); k++) //KfpDqSubdomains[linearSystems[0].ID].GetLength(0)=Mac
            {
                #region Set proper RHS 
                //var globalRHS = new Vector(model.TotalDOFs); //TODO: uncoomment if globalRHS is needed for solver
                foreach (ILinearSystem_v2 secSubdomain in linearSystems.Values)
                {

                    secSubdomain.RhsVector = Vector.CreateFromArray(KfpDqSubdomains[secSubdomain.Subdomain.ID][k], false);
                    //secSubdomain.RhsVector = Vector.CreateFromArray(KfpDqSubdomains[secSubdomain.Subdomain.ID][k], true); Wste sigoura na mhn peiraxthei to double[]


                    //mappings[seclinearSystems.Select((v, i) => new { System = v, Index = i }).First(x => x.System.ID == secSubdomain.ID).Index].SubdomainToGlobalVector(subdomainRHS.Data, globalRHS.Data);
                    //TODO: uncoomment if globalRHS is needed for solver
                }
                #endregion

                #region Solve
                solver.Solve();
                #endregion

                #region Copy solution in output vector
                foreach (ILinearSystem_v2 secSubdomain in linearSystems.Values)
                {
                    f2_vectorsSubdomains[secSubdomain.Subdomain.ID][k] = secSubdomain.Solution.CopyToArray();
                }
                #endregion
            }
            #endregion

            return f2_vectorsSubdomains;
        }

        public static Dictionary<int, double[][]> CalculateKpfKffinverseKfpDqSubdomains(Dictionary<int, double[][]> f2_vectorsSubdomains, Model_v2 model, IElementMatrixProvider elementProvider, IScaleTransitions_v2 scaleTransitions, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, double[][]> f3_vectorsSubdomains = new Dictionary<int, double[][]>();

            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                f3_vectorsSubdomains.Add(subdomain.ID, SubdomainCalculations_v2.CalculateKpfKffinverseKfpDq_v2(f2_vectorsSubdomains[subdomain.ID], subdomain, elementProvider, scaleTransitions, boundaryNodes));
            }
            return f3_vectorsSubdomains;
        }
        #endregion
    }
}
