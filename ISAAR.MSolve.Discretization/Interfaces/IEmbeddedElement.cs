using ISAAR.MSolve.Discretization.Embedding;
using System.Collections.Generic;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IEmbeddedElement
    {
        IList<EmbeddedNode> EmbeddedNodes { get; }
        Dictionary<DOFType, int> GetInternalNodalDOFs(IElement element, INode node);
        double[] GetLocalDOFValues(IElement hostElement, double[] hostDOFValues);
    }
}
