using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.PreProcessor.Elements
{
    static class PrintUtilities
    {
        public static void WriteToFile(double[,] array, string path)
        {
            var writer = new StreamWriter(path);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                for (int j = 0; j < array.GetLength(1); ++j)
                {
                    writer.Write(array[i, j]);
                    writer.Write(' ');
                }
                writer.WriteLine();
            }
            writer.Flush();
        }

        public static void WriteToFileVector(double[] array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                    writer2.Write(array[i]);
                    writer2.Write(' ');
                    writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }
            writer2.Flush();

        }

    }
}
