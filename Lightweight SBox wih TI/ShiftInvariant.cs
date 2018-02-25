using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using System.Threading;
using System.Collections.Concurrent;
namespace Lightweight_SBox_wih_TI
{
    class ShiftInvariant
    {
        public int size;//bit宽度
        public int degree;//布尔函数的次数
        //构造方法
        public ShiftInvariant(int s,int d)
        {
            size=s;
            degree=d;
        }
        //工具函数
        #region
        //小数量的阶乘函数
        int factorial(int x)
        {
            int result = 1;
            for (int i = 1; i <= x; i++)
                result *= i;
            return result;
        }
        //小数量的组合数
        int Comb(int n, int k)
        {
            int result = 1;
            for (int i = 0; i < k; i++)
                result = result * (n - i);
            return result / factorial(k);
        }
        //取bit操作
        int GetBit(long input, int pos)
        {
            return ((input >> pos) & (long)0x1) > 0 ? 1 : 0;
        }

        //取bit操作
        int GetBit(byte[] input, int pos)
        {
            int ind = pos / 8;
            int p = pos % 8;
            return (input[ind] >> p) & 0x01;
        }
        //置bit操作
        void SetBit(byte[] input, int pos)
        {
            int ind = pos / 8;
            int p = pos % 8;
            input[ind] = (byte)(input[ind] | (0x1 << p));
        }
        //置bit操作
        int SetBit(int input, int pos)
        {
            return (input| (0x1 << pos));
        }
        //循环右移1位
        void ShiftRight1(int varNum, int[] table)
        {
            int limit=(0x1<<varNum)-1;
            for(int i=0;i<table.Length;i++)
            {
                table[i] = (table[i] >> 1) | ((table[i] & 0x01) << (varNum - 1));
                table[i] = table[i] & limit;
            }
        }
        //由ANF计算真值表
        byte[] InvMoebiusTrans(int varNum, int[] ANF)
        {
            int[] ANF1 = new int[ANF.Length];
            for (int i = 0; i < (int)Math.Pow(2, varNum); i++)
                ANF1[i] = ANF[i];
            //分别定义small table size 和 small table position
            int sz, pos;

            for (int i = 0; i < varNum; i++)
            {
                sz = (int)Math.Pow(2, i);
                pos = 0;

                while (pos < (int)Math.Pow(2, varNum))
                {
                    for (int j = 0; j < sz; j++)
                    {
                        ANF1[pos + sz + j] = (ANF1[pos + j] ^ ANF1[pos + sz + j]);
                    }

                    pos = (pos + 2 * sz);
                }
            }
            int bytelen = (int)Math.Ceiling(Math.Pow(2, varNum) / 8);
            byte[] TruthTable=new byte[bytelen];
            for (int i = 0; i < bytelen; i++)
                TruthTable[i] = 0;
            for (int i = 0; i < (int)Math.Pow(2, varNum); i++)
            {
                if (ANF1[i] == 1)
                    SetBit(TruthTable, i);
            }
            return TruthTable;
        }

        //检查平衡性
        bool CheckBalanced(int varNum, byte[] TT)
        {
            int count = 0;
            for (int i = 0; i < ((int)0x1 << varNum); i++)
            {
                if (GetBit(TT, i) == 1)
                    count++;
            }
            if (count * 2 == ((int)0x1 << varNum))
                return true;
            else
                return false;
        }

        //检查置换性
        bool CheckPerm(int varNum, byte[] TT)
        {
            int len=0x1<<varNum;
            int[] input = new int[len];
            int[] output = new int[len];
            for (int i = 0; i < len; i++)
            {
                input[i] = i;
                output[i] = 0;
            }
            for(int j=0;j<varNum;j++)
            {
                for (int i = 0; i < len; i++)
                {
                    if (GetBit(TT, input[i]) == 1)
                        output[i] = SetBit(output[i], j);
                }
                //循环右移1位
                if (j != varNum - 1)
                    ShiftRight1(varNum, input);
            }
            
            //输出表排序
            Array.Sort(output);
            for (int i = 0; i < output.Length - 1; i++)
                if (output[i] == output[i + 1])
                    return false;
            return true;
        }

        //检查含有x0的项,且ANF[0]==0
        bool Checkx0(int[] ANF, int varNum)
        {
            if (ANF[0] == 1)
                return false;
            for (int i = 1; i < ANF.Length; i = i + 2)
                if (ANF[i] == 1)
                    return true;
            return false;
        }
        //计算ANF的次数
        int ANFdegree(int[] ANF,int v)
        {
            int maxd = 0;
            for (int i = 0; i < ANF.Length; i++)
            {
                if (ANF[i] == 0)
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
                if (d > maxd)
                    maxd = d;
            }
            return maxd;
        }
        #endregion


        //搜索满足条件的置换，将置换对应的真值表写入文件
        //统计平衡的个数和形成置换的个数
        public void SearchPermutation(string resultfile)
        {

            FileStream fs = new FileStream(resultfile, FileMode.Create);
            //StreamWriter sw = new StreamWriter(fs);
            BinaryWriter bw = new BinaryWriter(fs);
            //sw.WriteLine("degree={0},size={1}", degree, size);
            //搜索开始，首先估计总搜索空间
            int SearchLen=0;
            for (int i = 0; i <= degree; i++)
                SearchLen += Comb(size, i);
            long SearchSpace = (long)Math.Pow(2, SearchLen);
            //把每个项按HW排列，取所有<=2次的项,把项和序号之间的关系制表;之后根据RankTable可从序号找到相应的项
            int[] RankTable = new int[SearchLen];
            RankTable[0] = 0;
            HammingWeight hw = new HammingWeight(size);
            int j = 1;
            while (j < SearchLen)
            {
                //  hw.PrintState();
                RankTable[j] = hw.ReturnNum();
                hw.HwNext();
                j++;
            }
            //初始化循环中用到的数组
            int[] ANF = new int[(0x1<<size)];
            long num = 0;
            long count = 0;//总体满足条件的解
            long count_x0 = 0;//满足含有x0相关项
            long count_balanced = 0;//满足平衡的
            long count_perm = 0;//满足置换的
            //主循环
            for (num = 0; num < SearchSpace; num++)
            {
                //过程输出
                if ((num & 0xffffff) == 0)
                {
                    System.Console.WriteLine("num={0:x},x0F={1},BalancedF={2},PermF={3},Allf={4},percent={5}%", num,count_x0, count_balanced, count_perm, count,100*num/(double)SearchSpace);
                }
                //从num中获取ANF
                for (int i = 0; i < SearchLen; i++)
                {
                    if (GetBit(num, i) == 1)
                        ANF[RankTable[i]] = 1;
                    else
                        ANF[RankTable[i]] = 0;
                }
                //检查含有x0
                if (!Checkx0(ANF, size))
                    continue;
                count_x0++;
                //从ANF获得真值表
                byte[] TT= InvMoebiusTrans(size, ANF);

                //检查比特平衡性
                if (!CheckBalanced(size, TT))
                    continue;
                count_balanced++;
                //检查置换性
                if (!CheckPerm(size, TT))
                    continue;
                count_perm++;
                //密码性质检查


                //
                count++;
                //写出相关真值表
                //string hex = BitConverter.ToString(TT).Replace("-", string.Empty);
                //System.Console.WriteLine("TT={0}", hex);
               // sw.WriteLine("{0}", hex);
                bw.Write(TT);


            }
            //写出最终的统计结果
            System.Console.WriteLine("x0F={0},BalancedF={1},PermF={2},Allf={3}", count_x0, count_balanced, count_perm, count);
            //sw.WriteLine("x0F={0},BalancedF={1},PermF={2},Allf={3}", count_x0, count_balanced, count_perm, count);


            bw.Close();
            fs.Close();
        }

        public void WriteFlush(BinaryWriter bw, ConcurrentQueue<byte[]> cq,ref long count)
        {
            while (true)
            {
                if (cq.IsEmpty)
                    Thread.Sleep(10000);
                else
                {
                    while (cq.Count > 0)
                    {
                        byte[] TT = null;
                        cq.TryDequeue(out TT);
                        bw.Write(TT);
                        bw.Flush();
                        count = count + 1;
                    }
                }
            }
           
            //    bw.Flush();
        }
        //搜索满足TI条件的置换，将置换对应的真值表写入文件
        //统计平衡的个数和形成置换的个数
        public void SearchTIPermutation(string resultfile)
        {

            FileStream fs = new FileStream(resultfile, FileMode.Create);
            BinaryWriter bw = new BinaryWriter(fs);
            //搜索开始，首先估计总搜索空间
            int SearchLen = 0;
            for (int i = 0; i <= degree; i++)
                SearchLen += Comb(size, i);
            long SearchSpace = (long)Math.Pow(2, SearchLen);
            //把每个项按HW排列，取所有<=2次的项,把项和序号之间的关系制表;之后根据RankTable可从序号找到相应的项
            int[] RankTable = new int[SearchLen];
            RankTable[0] = 0;
            HammingWeight hw = new HammingWeight(size);
            int j = 1;
            while (j < SearchLen)
            {
                RankTable[j] = hw.ReturnNum();
                hw.HwNext();
                j++;
            }
            //初始化循环中用到的数组
           
            long count = 0;//总体满足条件的解
            long count_x0 = 0;//满足含有x0相关项
            long count_balanced = 0;//满足平衡的
            long count_perm = 0;//满足置换的
            long count_ti_balanced_bit = 0;//TI后bit满足平衡性的
            long count_ti_balanced_shares = 0;//TI后某bit的s个shares满足平衡性的
            long count_ti_perm = 0; //TI后满足置换性的
            ConcurrentQueue<byte[]> cq = new ConcurrentQueue<byte[]>();
            Thread td=new Thread(()=>{this.WriteFlush(bw,cq,ref count_ti_perm);});
            td.Start();
            //主循环
            Parallel.For(0, SearchSpace, new ParallelOptions { MaxDegreeOfParallelism = 10 }, num =>
            {
                int[] ANF = new int[(0x1 << size)];
                //过程输出
                if ((num & 0xfffff) == 0)
                {
                    System.Console.WriteLine("num={0:x},x0F={1},BalancedF={2},PermF={3},TI_BalancedF={6},TI_BalancedSharesF={7},TI_PermF={8},Allf={4},percent={5}%", num, count_x0, count_balanced, count_perm, count, 100 * num / (double)SearchSpace, count_ti_balanced_bit, count_ti_balanced_shares, count_ti_perm);
                }
                //从num中获取ANF
                for (int i = 0; i < SearchLen; i++)
                {
                    if (GetBit(num, i) == 1)
                        ANF[RankTable[i]] = 1;
                    else
                        ANF[RankTable[i]] = 0;
                }
                if (ANFdegree(ANF, size) < degree)
                    return;
                //检查含有x0
                if (!Checkx0(ANF, size))
                    return;
                Interlocked.Increment(ref count_x0);
                //从ANF获得真值表
                byte[] TT = InvMoebiusTrans(size, ANF);

                //检查比特平衡性
                if (!CheckBalanced(size, TT))
                    return;
                Interlocked.Increment(ref count_balanced);
                //检查置换性
                if (!CheckPerm(size, TT))
                    return;
                Interlocked.Increment(ref count_perm);
                //CheckPerm(size, TT);
                //检查TI
                SharedTransformation sf = new SharedTransformation(degree + 1, size, ANF);
                //检查TI后单bit的平衡性
                if (!sf.BitBalance())
                    return;
                //Interlocked.Increment(ref count_ti_balanced_bit);
                //检查TI后某个bit的S个share的平衡性
                if (!sf.SharedBitBalance())
                    return;
                Interlocked.Increment(ref count_ti_balanced_shares);
                if (!sf.Permutation())
                    return;

                //写出相关真值表
                cq.Enqueue(TT);
                //string hex = BitConverter.ToString(TT).Replace("-", string.Empty);
                //Interlocked.Increment(ref count_ti_perm);
                //Interlocked.Increment(ref count);
                //lock (bw)
                //{
                //    //System.Console.WriteLine("d={1},TT={0}", hex, ANFdegree(ANF, size));
                //    //sw.WriteLine("{0}", hex);
                //    bw.Write(TT);
                //    bw.Flush();
                //}

            });
            Thread.Sleep(10000);
            td.Abort();
            //写出最终的统计结果
            System.Console.WriteLine("x0F={0},BalancedF={1},PermF={2},TI_BalancedF={4},TI_BalancedSharesF={5},TI_PermF={6},Allf={3}", count_x0, count_balanced, count_perm, count_ti_perm, count_ti_balanced_bit, count_ti_balanced_shares, count_ti_perm);
            //sw.WriteLine("x0F={0},BalancedF={1},PermF={2},Allf={3}", count_x0, count_balanced, count_perm, count);
           
            bw.Close();
            fs.Close();
        }
    }
}
