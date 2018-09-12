using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IElement
    {
	    int ID { get; set; }
		IElementType IElementType { get; }  //TODO: The name should be just ElementType
        IList<INode> INodes { get; } //TODO: The name should be just Nodes
    }
}
