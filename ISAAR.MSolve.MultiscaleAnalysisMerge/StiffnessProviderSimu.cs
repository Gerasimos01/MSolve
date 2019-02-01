using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.MultiscaleAnalysisMerge
{
    public class StiffnessProviderSimu : IElementMatrixProvider
    {
        #region IElementMatrixProvider Members
        private SubdomainCalculationsSimultaneousObje host;

        public StiffnessProviderSimu(SubdomainCalculationsSimultaneousObje host)
        {
            this.host = host;
        }

        public IMatrix2D Matrix(IElement element)
        {
            var elementMatrix = element.IElementType.StiffnessMatrix(element);
            host.UpdateVectors(element, elementMatrix);
            return elementMatrix;
        }

        #endregion
    }
}






    
        //#region IElementMatrixProvider Members

        //public IMatrix2D Matrix(IElement element)
        //{
        //    return element.IElementType.StiffnessMatrix(element);
        //}

        //#endregion
    
