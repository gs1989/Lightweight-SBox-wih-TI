using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Lightweight_SBox_wih_TI
{
    //internal class
    //内部类：ANF表示转换为对象
    public class ANFterm: IComparable<ANFterm>
    {
        public int d;//degree
        public int base_index;//minimal bit index
        public int[] other_index;//all other relative index(index-base_index)
        public ANFterm(int[] indexes)
        {
            d = indexes.Length;
            //ascending?
            Array.Sort(indexes);
            base_index = indexes[0];
            if (d > 1)
            {
                other_index = new int[d - 1];
                for (int i = 0; i < d - 1; i++)
                    other_index[i] = indexes[i + 1] - indexes[0];
            }
        }
        public int CompareTo(ANFterm t1)
        {
            if (t1.d != this.d)
                return this.d - t1.d;
            else
            {
                if(t1.d==1)
                    return this.base_index - t1.base_index;
                if (this.other_index[0] != t1.other_index[0])
                    return this.other_index[0] - t1.other_index[0];
                else
                {
                    return this.base_index - t1.base_index;
                }
            }
        }
    }
}
