using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Loading.Interfaces;

namespace ISAAR.MSolve.IGA.Loading.NodalLoads
{
    public class NodalLoad :Load, ILoad
    {
        public NodalLoad(INode node, IDofType dofType, double amount)
        {
            Node = node;
            DOF = dofType;
            Amount = amount;
        }

        public Table<INode, IDofType, double> CalculateLoad()
        {
            var loadTable = new Table<INode, IDofType, double>();
            loadTable.TryAdd(Node, DOF, Amount);
            return loadTable;
        }
    }
}
