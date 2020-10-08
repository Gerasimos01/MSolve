using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.IGA.Loading.Interfaces
{
    public interface ILoad
    {
        Table<INode, IDofType, double> CalculateLoad();
    }
}