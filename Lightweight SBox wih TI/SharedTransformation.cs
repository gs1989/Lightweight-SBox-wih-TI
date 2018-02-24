using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Lightweight_SBox_wih_TI
{
    class SharedTransformation
    {
        public int s;//share数量
        public int v;//变量个数
        public int[] TT1;//bit 0中share0对应的真值表

        public SharedTransformation(int si, int vi, int[] ANF)//输入share前的ANF
        {
            s = si;
            v = vi;
            TI_wrapper tw=new TI_wrapper(vi, si,ANF);
            tw.Compute_From_ANF();
            TT1 = tw.Get_TruthTable_F1();
        }

        //检验是否单bit是否平衡
        public bool BitBalance()
        {
            int sum=0;
            for (int i = 0; i < TT1.Length; i++)
                sum += TT1[i];
            if ((2 * sum) == TT1.Length)
                return true;
            else
                return false;
        }
        //由完整s*v bit信息计算F1的输出
        //!!!!!!这里的字节序可能不对，只是用来检验平衡性和置换性，因此不用修正
        //不能用来输出真值表!!!
        public int ComputeF1(int x)
        {
            int t = 0;
            int mask=(0x1<<(s-1))-1;
            for (int i = 0; i < v; i++)
            {
                t =t| (x & mask);
                if (i == v - 1)
                    break;
                x = x >> s;
                t = t << (s-1);
            }
            return TT1[t];

        }
        //Shares向右循环移位1
        public int ShiftShares(int x)
        {
            int mask = (0x1 << s)-1;
            int[] xb = new int[v];
            //split to v variables
            for(int i=0;i<v;i++)
            {
                xb[i] = x & mask;
                x = x >> s;
            }
            //Shift Right 1
            for (int i = 0; i < v; i++)
            {
                xb[i] = (xb[i] >> 1) | ((xb[i] & 0x01) << (s - 1));
            }
            //paste back to one
            int y = 0;
            for (int i = v - 1; i > -1; i--)
            {
                y = y | xb[i];
                if(i!=0)
                   y = y << s;
            }
            return y;
        }
        //原始输入bit向右循环移位1
        public int ShiftBit(int x)
        {
            int mask = (0x1 << s) - 1;
            int[] xb = new int[v];
            int[] yb = new int[v];
            //split to v variables
            for (int i = 0; i < v; i++)
            {
                xb[i] = x & mask;
                x = x >> s;
            }
            //Shift Right 1
            for (int i = 0; i < v; i++)
            {
                yb[i] = xb[(i + 1) % v];
            }
            //paste back to one
            int y = 0;
            for (int i = v - 1; i > -1; i--)
            {
                y = y | yb[i];
                if (i != 0)
                    y = y << s;
            }
            return y;
        }
        //检验bit 0的s个share中是否平衡
        public bool SharedBitBalance()
        {
            int limit = 0x1 << (s * v);
            int[] count=new int[0x1<<s];
            int result = 0;
            int xt = 0;
            for (int x = 0; x < limit; x++)
            {
                result = 0;
                xt = x;
                for (int i = 0; i < s; i++)
                {
                    //计算当前F1的输出
                    result |= ComputeF1(xt);
                    if (i != s - 1)
                    {
                        //移动x
                        xt = ShiftShares(xt);
                        //移动result
                        result = result << 1;
                    }
                }
                count[result]++;
            }
            for (int i = 1; i < count.Length; i++)
                if (count[i] != count[0])
                    return false;
            return true;
        }

        //检验share后的大置换是否仍为置换
        public bool Permutation()
        {
            int limit = 0x1 << (s * v);
            int[] Bit1Table = new int[limit];
            int result = 0;
            int xt = 0;
            //Generate bit 1 Table
            for (int x = 0; x < limit; x++)
            {
                result = 0;
                xt = x;
                for (int i = 0; i < s; i++)
                {
                    //计算当前F1的输出
                    result |= ComputeF1(xt);
                    if (i != s - 1)
                    {
                        //移动x
                        xt = ShiftShares(xt);
                        //移动result
                        result = result << 1;
                    }
                }
                Bit1Table[x] = result;
            }
            //根据Bit1Table算出剩余的置换
            bool[] Used=new bool[limit];
            for (int x = 0; x < limit; x++)
            {
                result = 0;
                xt = x;
                for (int i = 0; i < v; i++)
                {
                    //计算当前F1的输出
                    result |= Bit1Table[xt];
                    if (i != v - 1)
                    {
                        //移动x
                        xt = ShiftBit(xt);
                        //移动result
                        result = result << s;
                    }
                }
                if (Used[result])
                    return false;
                else
                    Used[result] = true;
            }
            return true;
        }
    }
}
