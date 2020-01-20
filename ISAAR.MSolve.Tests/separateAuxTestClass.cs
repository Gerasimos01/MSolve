using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;
using System;
using ISAAR.MSolve.LinearAlgebra.Commons;
using System.Linq;





namespace ISAAR.MSolve.SamplesConsole
{
    public static class separateAuxTestClass
    {
        [Fact]
        public static void Check01()
        {
            var dictionary1 = new Dictionary<int, double[]>() { { 1, new double[2] }, { 2, new double[2] }, { 3, new double[2] }, };
            var dictionary2 = new Dictionary<int, double[]>() { { 1, new double[2] }, { 2, new double[2] }, };

            var diferrentKeys = dictionary1.Keys.Where(x => !dictionary2.ContainsKey(x)).ToList();
            
            

        }
    }
}
