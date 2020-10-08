using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.IGA.Elements.Boundary;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.IGA.Postprocessing
{
    public class ParaviewIsogeometricSolids
    {
        private readonly Model _model;
        private readonly IVectorView _solution;
        private readonly string _filename;

        public ParaviewIsogeometricSolids(Model model, IVectorView solution, string filename)
        {
            _model = model;
            _solution = solution;
            _filename = filename;
        }

        public void CreateParaviewFile()
        {
            throw new NotImplementedException();
        }

        private int[,] CreateConnectivity()
        {
			throw new NotImplementedException();
        }

        private double[,] CalculateProjectiveControlPoints()
        {
            var projectiveCPs = new double[_model.PatchesDictionary[0].ControlPoints.Count, 4];
            foreach (var controlPoint in _model.PatchesDictionary[0].ControlPoints)
            {
                var weight = controlPoint.WeightFactor;
                projectiveCPs[controlPoint.ID, 0] = controlPoint.X * weight;
                projectiveCPs[controlPoint.ID, 1] = controlPoint.Y * weight;
                projectiveCPs[controlPoint.ID, 2] = controlPoint.Z * weight;
                projectiveCPs[controlPoint.ID, 3] = weight;
            }

            return projectiveCPs;
        }

        private void WriteParaviewFile(double[,] nodeCoordinates, int[,] elementConnectivity, double[,] displacements)
		{
			
		}
    }
}
