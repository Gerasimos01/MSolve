using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.MultiscaleAnalysisMerge
{
    public class StiffnessProviderSimu_v2 : IElementMatrixProvider
    {
        #region IElementMatrixProvider Members
        private SubdomainCalculationsSimultaneousObje_v2 host;

        public StiffnessProviderSimu_v2(SubdomainCalculationsSimultaneousObje_v2 host)
        {
            this.host = host;
        }

        public IMatrix2D Matrix(IElement element) //TODOGer IMatrix2D will be changed to Matrix etc.
        {
            var elementMatrix = element.IElementType.StiffnessMatrix(element);
            host.UpdateVectors_v2(element, elementMatrix);
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
    
