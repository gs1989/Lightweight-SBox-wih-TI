using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Lightweight_SBox_wih_TI
{
    class HammingWeight
    {
        public bool[] M;
        public int dim;
        public int pstart;
        public int pend;
        public int Hw;
        public HammingWeight(int d)
        {
            int i;
            M = new bool[d];
            dim = d;
            pstart = 0;
            pend = 0;
            Hw = 1;
            for (i = 1; i < dim; i++)
                M[i] = false;
            M[0] = true;
        }
        public HammingWeight(int d, int m)
        {
            int i;
            M = new bool[d];
            dim = d;
            pstart = 0;
            pend = m - 1;
            Hw = 1;
            for (i = m; i < dim; i++)
                M[i] = false;
            for (i = 0; i < m; i++)
                M[i] = true;
        }
        //更新向量
        public void HwNext()
        {
            int i;

            //更新之后，Hw不必变化
            if (pstart != dim - Hw)
            {
                if (pstart == pend)//考虑特殊情况
                {
                    M[pstart + 1] = true;

                    for (i = 0; i < pstart + 1; i++)
                    {
                        M[i] = false;
                    }

                    //更新参数
                    pend = pstart + 1;
                    pstart = pstart + 1;

                    while (((pend + 1) < dim) && (M[pend + 1] != false))
                    {
                        pend = pend + 1;
                    }

                }

                else  //pstart != pend
                {
                    M[pend + 1] = true;

                    for (i = 0; i < pend - pstart; i++)
                    {
                        M[i] = true;
                    }

                    for (i = pend - pstart; i < pend + 1; i++)
                    {
                        M[i] = false;
                    }

                    //更新pstart和pend
                    pend = (pend - pstart - 1);
                    pstart = 0;
                }
            }

            //出现00...0011...11的情况
            else
            {
                if (Hw == dim)
                {
                    Hw++;
                    return;
                }
                for (i = 0; i < Hw + 1; i++)
                {
                    M[i] = true;
                }

                for (i = Hw + 1; i < dim; i++)
                {
                    M[i] = false;
                }

                //更新参数
                pstart = 0;
                pend = Hw;
                Hw = Hw + 1;

            }
        }

        public void PrintState()
        {
            System.Console.Write("M=");
            for (int i = 0; i < dim; i++)
                System.Console.Write("{0}\t", (M[i]) ? 1 : 0);
            System.Console.Write("Hw={0}\n", Hw);
        }

        //返回对应序号数组
        public int[] ReturnNo()
        {
            int[] No = new int[Hw];
            int i, j = 0;
            for (i = 0; i < dim; i++)
            {
                if (M[i])
                {
                    No[j] = i;
                    j++;
                    if (j == Hw)
                        break;
                }
            }
            return No;
        }
        //返回对应序号数
        public int ReturnNum()
        {
            int Num =0;
            int i = 0;
            for (i = dim-1; i >-1; i--)
            {
                if (M[i])
                {
                    Num++;
                }
                if(i!=0)
                Num = Num << 1;
            }
            return Num;
        }
    }
}
