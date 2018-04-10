using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
namespace Lightweight_SBox_wih_TI
{
    class TI_wrapper
    {
        public int v;//变元个数
        public int[] Ori_ANF;//原始ANF
        public int s;//share个数
        public int[] TI_Form;//TI表达式，定义为一个长为v的s+1进制数，为0表示无此变元，为其他表示此变元的下标
        public TI_wrapper(int vi, int si, int[] ANF)
        {
            v = vi;
            s = si;
            Ori_ANF = new int[ANF.Length];
            for (int i = 0; i < ANF.Length; i++)
                Ori_ANF[i] = ANF[i];
            TI_Form = new int[(int)Math.Pow(s + 1, v)];
        }
        public int GetBit(int x, int no)
        {
            int mask = 0x01;
            mask = mask << no;
            if ((x & mask) > 0)
                return 1;
            else
                return 0;
        }
        public int SetBit(int x, int no)
        {
            int mask = 0x01;
            mask = mask << no;
            return (x | mask);
        }
        //将一个整数分解为d长的s进制数
        public int[] Trans(int num, int s, int d)
        {
            int[] t = new int[d];
            int temp = num;
            int i = 0;
            while (temp != 0)
            {
                t[i] = temp % (s);
                i++;
                temp = temp / (s);
            }
            return t;
        }
        public int Compute_TruthValue(int[][] vshared)
        {
            int result = 0;
            for (int i = 0; i < TI_Form.Length; i++)
            {
                if (TI_Form[i] == 1)
                {
                    int[] temp = Trans(i, s + 1, v);
                    int tempresult = -1;
                    for (int j = 0; j < v; j++)
                        if (temp[j] != 0)//等于0的说明缺项，不用管
                        {
                            if (tempresult == -1)
                                tempresult = 1;
                            tempresult *= vshared[j][temp[j] - 1];
                        }
                    if (tempresult == -1)//全都缺项，说明为常数项
                        tempresult = 1;
                    result = result ^ tempresult;
                }
            }
            return result;
        }
        //真值表转ANF
        public int[] TruthTable2ANF(int[] tt)
        {
            int n = (int)Math.Log(tt.Length, 2);
            int[] ANF = new int[tt.Length];
            for (int i = 0; i < tt.Length; i++)
                ANF[i] = tt[i];
            for (int i = 0; i < n; i++)
            {
                int Sz = (int)Math.Pow(2, i);
                int pos = 0;
                while (pos < Math.Pow(2, n))
                {
                    for (int j = 0; j < Sz; j++)
                        ANF[pos + Sz + j] = ANF[pos + j] ^ ANF[pos + Sz + j];
                    pos = pos + 2 * Sz;
                }
            }
            return ANF;
        }
        //s*v<32!!!!!!!
        //每个变元赋值，内部再分为2^s次种取值，分别计算TI_form的值，要求这2^s中取值均是相同的
        //返回整体对非share值的真值表
        public void Print_RealTruthTable()
        {
            int len = v * s;
            int max = (int)Math.Pow(2, len);
            int[][] var_shared = new int[v][];
            int[] var = new int[v];
            int[] TruthTable = new int[(int)Math.Pow(2, v)];
            for (int i = 0; i < v; i++)
            {
                var_shared[i] = new int[s];
            }
            for (int i = 0; i < TruthTable.Length; i++)
                TruthTable[i] = -1;
            System.Console.Write("真值表:\t");
            for (int input = 0; input < max; input++)//每个输入取值
            {
                int outputindex = 0;
                for (int i = 0; i < v; i++)//取v个变元
                {
                    var[i] = 0;
                    for (int j = 0; j < s; j++)
                    {
                        var_shared[i][j] = GetBit(input, i * s + j);
                        var[i] ^= var_shared[i][j];
                    }
                    if (var[i] == 1)
                        outputindex = SetBit(outputindex, i);//v1在低位
                }
                int newv = Compute_TruthValue(var_shared);
                if (TruthTable[outputindex] == -1)
                    TruthTable[outputindex] = newv;
                else
                    if (TruthTable[outputindex] != newv)
                    {
                        System.Console.WriteLine("不同share分配返回结果不同!");
                    }

            }
            for (int input = 0; input < TruthTable.Length; input++)
                System.Console.Write("{0}\t", TruthTable[input]);
            System.Console.WriteLine();
            int[] ANF = TruthTable2ANF(TruthTable);
            System.Console.Write("ANF:\t");
            for (int input = 0; input < ANF.Length; input++)
                System.Console.Write("{0}\t", ANF[input]);
            System.Console.WriteLine();
            for (int i = 0; i < TruthTable.Length; i++)
                if (ANF[i] != Ori_ANF[i])
                    System.Console.WriteLine("出错！与原ANF不符...");
        }

        public void Print_TruthTable_F1(String filename, int num)
        {
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            int len = v * (s - 1);
            int max = (int)Math.Pow(2, len);
            int[][] var_shared = new int[v][];
            int[] var = new int[v];
            int[] TruthTable = new int[(int)Math.Pow(2, v * (s - 1))];
            for (int i = 0; i < v; i++)
            {
                var_shared[i] = new int[s];
            }
            for (int i = 0; i < TruthTable.Length; i++)
                TruthTable[i] = -1;
            for (int input = 0; input < max; input++)
            {
                for (int i = 0; i < v; i++)
                {
                    for (int j = 1; j < s; j++)
                    {
                        var_shared[i][j] = GetBit(input, i * (s - 1) + j - 1);
                    }
                }
                int newv = Compute_TruthValue(var_shared);
                TruthTable[input] = newv;
            }
            sw.WriteLine("module F_I{0}O{1}_D{2}_{3}(in,out);", v, 1, s - 1, num);
            sw.WriteLine("  input[{0}:0] in;", len - 1);
            sw.WriteLine("  output out;");
            sw.WriteLine("  reg out;");
            sw.WriteLine("  always@(in)");
            sw.WriteLine("  begin");
            sw.WriteLine("    case(in)");
            for (int i = 0; i < max; i++)
                sw.WriteLine("   {0}: out<={1};", i, TruthTable[i]);
            sw.WriteLine(" endcase");
            sw.WriteLine("end");
            sw.WriteLine("endmodule");
            sw.Close();
            fs.Close();
        }
        //返回真值表，以bit分组，每个组里为s个shares
        public int[] Get_TruthTable_F1()
        {
            int len = v * (s - 1);
            int max = (int)Math.Pow(2, len);
            int[][] var_shared = new int[v][];
            int[] var = new int[v];
            int[] TruthTable = new int[(int)Math.Pow(2, v * (s - 1))];
            for (int i = 0; i < v; i++)
            {
                var_shared[i] = new int[s];
            }
            for (int i = 0; i < TruthTable.Length; i++)
                TruthTable[i] = -1;
            for (int input = 0; input < max; input++)
            {
                //按bit排序
                for (int i = 0; i < v; i++)
                {
                    for (int j = 1; j < s; j++)
                    {
                        var_shared[i][j] = GetBit(input, i * (s - 1) + j - 1);
                    }
                }
                
                int newv = Compute_TruthValue(var_shared);
                TruthTable[input] = newv;
            }
            return TruthTable;
        }
        //返回真值表，以shares分组，每个分组为v个比特
        public int[] Get_TruthTable_F1_SharesGroup()
        {
            int len = v * (s - 1);
            int max = (int)Math.Pow(2, len);
            int[][] var_shared = new int[v][];
            int[] var = new int[v];
            int[] TruthTable = new int[(int)Math.Pow(2, v * (s - 1))];
            for (int i = 0; i < v; i++)
            {
                var_shared[i] = new int[s];
            }
            for (int i = 0; i < TruthTable.Length; i++)
                TruthTable[i] = -1;
            for (int input = 0; input < max; input++)
            {
                //按bit排序
                for (int i = 0; i < v; i++)
                {
                    for (int j = 1; j < s; j++)
                    {
                        var_shared[i][j] = GetBit(input, (j-1) * v +i);//第i bit第j个shares，位置应为(j-1)*v+i
                    }
                }

                int newv = Compute_TruthValue(var_shared);
                TruthTable[input] = newv;
            }
            return TruthTable;
        }
        public void Print_TruthTable_F1(String filename, long num)
        {
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            int len = v * (s - 1);
            int max = (int)Math.Pow(2, len);
            int[][] var_shared = new int[v][];
            int[] var = new int[v];
            int[] TruthTable = new int[(int)Math.Pow(2, v * (s - 1))];
            for (int i = 0; i < v; i++)
            {
                var_shared[i] = new int[s];
            }
            for (int i = 0; i < TruthTable.Length; i++)
                TruthTable[i] = -1;
            for (int input = 0; input < max; input++)
            {
                for (int i = 0; i < v; i++)
                {
                    for (int j = 1; j < s; j++)
                    {
                        var_shared[i][j] = GetBit(input, i * (s - 1) + j - 1);
                    }
                }
                int newv = Compute_TruthValue(var_shared);
                TruthTable[input] = newv;
            }
            sw.WriteLine("module F_I{0}O{1}_D{2}_{3:x}(in,out);", v, 1, s - 1, num);
            sw.WriteLine("  input[{0}:0] in;", len - 1);
            sw.WriteLine("  output out;");
            sw.WriteLine("  reg out;");
            sw.WriteLine("  always@(in)");
            sw.WriteLine("  begin");
            sw.WriteLine("    case(in)");
            for (int i = 0; i < max; i++)
                sw.WriteLine("   {0}: out<={1};", i, TruthTable[i]);
            sw.WriteLine(" endcase");
            sw.WriteLine("end");
            sw.WriteLine("endmodule");
            sw.Close();
            fs.Close();
        }
        //从ANF开始计算TI
        public void Compute_From_ANF()
        {
            for (int i = 0; i < Ori_ANF.Length; i++)
            {
                if (Ori_ANF[i] == 0)
                    continue;
                //对原始ANF中的每一项
                int[] con = new int[v];
                int d = 0;
                for (int j = 0; j < v; j++)
                {
                    if (GetBit(i, j) == 1)//修改
                    {
                        con[j] = 1;
                        d++;
                    }
                }
                //d表示次数，con表示其中含有的bit
                if (d == 0)
                {
                    TI_Form[0] ^= 1;
                    continue;
                }
                TI_Searcher ts = new TI_Searcher(d, s);
                ts.ComputeCoef();
                //2018.3.26 Ensure y0 contains x0,x1, do not have x2
                //ts.Shift1();


                int[] tempTI = ts.GetCorrectTerm(con, v);



                //验证用，加和所有share---------------
                //ts.Shift1();
                //for (int k = 1; k < s; k++)
                //{
                //    int[] temp = ts.GetCorrectTerm(con, v);
                //    for (int j = 0; j < tempTI.Length; j++)
                //    {
                //        tempTI[j] ^= temp[j];
                //    }
                //    ts.Shift1();
                //}
                //实际使用中应该去掉这部分---------------

                for (int j = 0; j < tempTI.Length; j++)
                {
                    TI_Form[j] ^= tempTI[j];
                }   


            }
            //Print_RealTruthTable();
            //System.Console.WriteLine("最终的TI和形式为");
            //int sum = 0;
            //for (int i = 0; i < TI_Form.Length; i++)
            //{
            //    System.Console.Write("{0}\t", TI_Form[i]);
            //    sum += TI_Form[i];
            //}
            //System.Console.WriteLine("\n共有{0}项", sum);
        }

        public void Test()
        {
            TI_Searcher ts = new TI_Searcher(2, s);
            ts.ComputeCoef();
            int[] con = { 1, 1 };
            for (int k = 0; k < s; k++)
            {
                int[] tempTI = ts.GetCorrectTerm(con, 2);
                for (int i = 0; i < tempTI.Length; i++)
                {
                    if (TI_Form[i] == 1 && tempTI[i] == 1)
                        System.Console.WriteLine("Error!");
                    TI_Form[i] ^= tempTI[i];
                }
                ts.Shift1();
            }
            System.Console.WriteLine("最终a*b对应的TI和形式为");
            int sum = 0;
            for (int i = 0; i < TI_Form.Length; i++)
            {
                System.Console.Write("{0}\t", TI_Form[i]);
                sum += TI_Form[i];
            }
            System.Console.WriteLine("\n共有{0}项", sum);

            Print_RealTruthTable();
        }

    }
}
