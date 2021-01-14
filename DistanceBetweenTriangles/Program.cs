using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Media3D;

using System.Diagnostics;

using Geometry3D;

namespace DistanceBetweenTriangles
{
    class Program
    {
        static void Main(string[] args)
        {
            //sizeof(Triangle3D);
            var ps = System.Runtime.InteropServices.Marshal.SizeOf(typeof(Point3D));
            var ls = System.Runtime.InteropServices.Marshal.SizeOf(typeof(LineSegment3D));
            var ts = System.Runtime.InteropServices.Marshal.SizeOf(typeof(Triangle3D));
        }
    }
}
