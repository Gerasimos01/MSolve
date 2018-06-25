using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;


namespace ISAAR.MSolve.MultiscaleAnalysis.Interfaces
{
    public interface IRVEbuilder
    {
        Tuple<Model, Dictionary<int, Node>,double> GetModelAndBoundaryNodes();
    }
}
