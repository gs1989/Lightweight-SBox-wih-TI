using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using System.Threading;
namespace Lightweight_SBox_wih_TI
{
    class Searcher8
    {
        public int size;//输入长度
        public int round;//迭代轮数
        public int len;//表长
        public int shift;//2bit输出之间输入的移位次数
        FileStream fsScript;
        StreamWriter swScript;
        int[] RankTable;//2次以下项的排序列表，存储ANF项的标号

        public Searcher8(int s,int r,int sh)
        {
            size = s;
            round = r;
            len = (int)Math.Pow(2, size);
            shift = sh;
            //生成对应排序列表
            int nsize = size / 4 * 3;
            int m=1+nsize+nsize*(nsize-1)/2;
            RankTable=new int[m];
            RankTable[0] = 0;
            HammingWeight hw = new HammingWeight(nsize);
            int i = 1;
            while (i < m)
            {
              //  hw.PrintState();
                RankTable[i] = hw.ReturnNum();
                hw.HwNext();
                i++;
            }
        }
        /// <summary>
        /// 取比特函数：取input的第pos位，其中第0位处于最右边，即LSB
        /// </summary>
        /// <param name="input"></param>
        /// <param name="pos"></param>
        /// <returns></returns>
        public int GetBit(int input, int pos)
        {
            return ((input >> pos) & 0x1);
        }
        public int GetBit(byte[] input, int pos)
        {
            return ((input[pos/8] >> (pos%8)) & 0x1);
        }
        public int GetBit(long input, int pos)
        {
            return ((input >> pos) & (long)0x1)>0?1:0;
        }
        //置bit操作
        int SetBit(int input, int pos)
        {
            return (input | (0x1 << pos));
        }
        /// <summary>
        /// 置比特函数：取input的第pos位，其中第0位处于最右边，即LSB
        /// </summary>
        /// <param name="input"></param>
        /// <param name="pos"></param>
        /// <returns></returns>
        public long SetBit(long input, int pos)
        {
            return (input | ((long)0x1 << pos));
        }

        /// <summary>
        /// 绝对值函数，返回input的绝对值
        /// </summary>
        /// <param name="input"></param>
        /// <returns></returns>
        public int Abs(int input)
        {
            if (input < 0)
            {
                return ((-1) * input);
            }
            else
            {
                return input;
            }
        }

        /// <summary>
        /// 幂次函数：计算baseNum^expNum
        /// </summary>
        /// <param name="baseNum"></param>
        /// <param name="expNum"></param>
        /// <returns></returns>
        public long Power(int baseNum, int expNum)
        {
            if (expNum == 0)
            {
                return 1;
            }
            else
            {
                return (long)(baseNum * Power(baseNum, expNum - 1));
            }
        }
        /// <summary>
        /// 计算S盒的差分均匀度，sbox是其真值表表示，size是其比特级规模
        /// </summary>
        /// <param name="sbox"></param>
        /// <param name="size"></param>
        /// <returns></returns>
        ///         /// <summary>
        /// 对于大小为（2^size）的table，返回其中的最大值
        /// </summary>
        /// <param name="table"></param>
        /// <param name="size"></param>
        /// <returns></returns>
        public int MaxTable(int[] table)
        {
            int maxValue = 0;

            for (int i = 0; i < len; i++)
            {
                if (table[i] > maxValue)
                {
                    maxValue = table[i];
                }
            }
            return maxValue;
        }
        public int DiffUniform(int[] sbox, int size)
        {
            int maxDiff;
            len = (int)Math.Pow(2, size);
            //差分分布表中的一行，即输入差分不变，输出差分根据输入值而改变
            int[] consInVarOut = new int[len];

            //差分分布表一共（2^size）行，记录每一行的最大值
            int[] multiRows = new int[len];

            for (int inDiff = 1; inDiff < len; inDiff++)
            {
                //对于某一输入差分inDiff，清空consInVarOut
                for (int i = 0; i < len; i++)
                {
                    consInVarOut[i] = 0;
                }

                //遍历输入值input，求解输出差分out
                for (int input = 0; input < len; input++)
                {
                    int outDiff;
                    outDiff = (sbox[input] ^ sbox[input ^ inDiff]);
                    consInVarOut[outDiff]++;
                }

                //做完DDT的一行之后，取出这行的最大值
                multiRows[inDiff] = MaxTable(consInVarOut);
            }

            maxDiff = MaxTable(multiRows);

            //释放consInVarOut
            //释放multiRows

            return maxDiff;
        }
        /// <summary>
        /// 内积函数：输入为2个int型数据，bitLength为内积的长度，从0计数
        /// </summary>
        /// <param name="bitLength"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public int InnerMulti(int bitLength, int a, int b)
        {
            int sum = 0;

            for (int i = 0; i < bitLength; i++)
            {
                sum = (sum ^ (GetBit(a, i) & GetBit(b, i)));
            }
            return sum;
        }
        /// <summary>
        /// 计算非线性度：sbox为S盒的真值表表示，size为其比特级规模
        /// </summary>
        /// <param name="sbox"></param>
        /// <param name="size"></param>
        /// <returns></returns>
        public int Nonlinear(int[] sbox,int size)
        {
            //一系列Walsh谱值,由于遍历u，v，所以一共大约2^(2*size)个Walsh谱绝对值
            int[] walsh = new int[Power(2, 2 * size)];
            int counter = 0;

            //输入掩码、输出掩码和输入值均为size长的0-1向量，直接用int型来表示
            for (int u = 0; u < len; u++)
            {
                for (int v = 1; v < len; v++)
                {
                    //记录偶数项和奇数项的总量
                    int oddNum = 0;
                    int evenNum = 0;

                    for (int input = 0; input < len; input++)
                    {
                        if (((InnerMulti(size, v, sbox[input]) ^ InnerMulti(size, u, input)) % 2) == 0)
                        {
                            evenNum++;
                        }
                        else
                        {
                            oddNum++;
                        }
                    }

                    //针对特定的(u,v)，遍历输入input后，得到对应的Walsh谱值(绝对值)
                    walsh[counter++] = Abs(evenNum - oddNum);
                }
            }

            //计算最大的walsh谱绝对值
            int maxValue = 0;
            for (int i = 0; i < Power(2, 2 * size); i++)
            {
                if (walsh[i] > maxValue)
                {
                    maxValue = walsh[i];
                }
            }

            //释放walsh
            return (int)(Power(2, size - 1) - maxValue / 2);
        }

        /// <summary>
        /// 由布尔函数的真值表计算其ANF
        public void MoebiusTrans(int varNum, int func, int[] ANF)
        {
            for (int i = 0; i < Power(2, varNum); i++)
            {
                ANF[i] = GetBit(func, i);
            }

            //分别定义small table size 和 small table position
            int sz, pos;

            for (int i = 0; i < varNum; i++)
            {
                sz = (int)Power(2, i);
                pos = 0;

                while (pos < Power(2, varNum))
                {
                    for (int j = 0; j < sz; j++)
                    {
                        ANF[pos + sz + j] = (ANF[pos + j] ^ ANF[pos + sz + j]);
                    }

                    pos = (pos + 2 * sz);
                }
            }
        }

        /// <summary>
        /// 由ANF计算其真值表
        public long InvMoebiusTrans(int varNum, int[] ANF)
        {
            int[] ANF1 = new int[ANF.Length];
            for (int i = 0; i < Power(2, varNum); i++)
                ANF1[i] = ANF[i];
            //分别定义small table size 和 small table position
            int sz, pos;

            for (int i = 0; i < varNum; i++)
            {
                sz = (int)Power(2, i);
                pos = 0;

                while (pos < Power(2, varNum))
                {
                    for (int j = 0; j < sz; j++)
                    {
                        ANF1[pos + sz + j] = (ANF1[pos + j] ^ ANF1[pos + sz + j]);
                    }

                    pos = (pos + 2 * sz);
                }
            }
            long func = 0;
            for (int i = 0; i < Power(2, varNum); i++)
            {
                if (ANF1[i] == 1)
                    func = SetBit(func, i);
            }
            return func;
        }

        //4分支，每个分支2比特
        public void OneRoundTrans(int[] table, long func)
        {
            int newX0, newX1, newX2, newX3;

            for (int i = 0; i < len; i++)
            {
                int Fin = table[i] >> 2;
                int Fout1 = GetBit(func, Fin);
                int mask = (0x01 << (shift)) - 1;
                int Fout2 = GetBit(func, ((Fin>>(shift))|((Fin&mask)<<(6-shift))));
                int Fout = (Fout2 << 1) | Fout1;
                int x0 = (GetBit(table[i], 1) << 1) | GetBit(table[i], 0);
                int x1 = (GetBit(table[i], 3) << 1) | GetBit(table[i], 2);
                int x2 = (GetBit(table[i], 5) << 1) | GetBit(table[i], 4);
                int x3 = (GetBit(table[i], 7) << 1) | GetBit(table[i], 6);
                newX3 = Fout ^ x0;
                newX2 = x3;
                newX1 = x2;
                newX0 = x1;
                table[i] = ((newX3 << 6) | (newX2 << 4) | (newX1 << 2) | newX0);
            }
        }

        //5-3卡方函数+拉线P+固定常数C
        public void OneRoundTrans_X2_BitP_FixC(int[] table, int[] Ptable,int C)
        {
            int xl,xr;//xl为5bit,xr为3bit
            int xl1, xl2, xr1, xr2;//循环移位结果
            int y,yl, yr;
            int pout;
            for (int i = 0; i < len; i++)
            {
                //拆分
                xl = (table[i]) >> 3;
                xr = table[i] & 0x07;
                //卡方函数
                xl1 = ((xl << 1) | (xl >> 4))&0x1f;
                xl2 = ((xl << 2) | (xl >> 3))&0x1f;
                xr1 = ((xr << 1) | (xr >> 2))&0x07;
                xr2 = ((xr << 2) | (xr >> 1))&0x07;

                yl = xl ^ ((xl1 ^ 0x1f) & xl2);
                yr = xr ^ ((xr1 ^ 0x07) & xr2);
                //合并
                y = (yl << 3) | yr;
                //bitP
                pout = 0;
                for (int j = 0; j < 8; j++)
                    if (GetBit(y, Ptable[j]) == 1)
                        pout=SetBit(pout, j);
                //加常数
                table[i] = pout^C;
            }
        }
        //两个4bit ShiftInvariant置换+8bit SI P
        public void OneRoundTrans_ShiftInvariant_SPN(int[] table, int[] Ptable,int[] Stable1,int[] Stable2)
        {
            int xl, xr;//xl为5bit,xr为3bit
            int Sout = 0;
            for (int i = 0; i < len; i++)
            {
                //拆分
                xl = (table[i]) >> 4;
                xr = table[i] & 0x0f;
                //S
                Sout = (Stable1[xl] << 4) | Stable2[xr];
                //P
                table[i] = Ptable[Sout];
            }
        }
        //5bit+3bit ShiftInvariant置换+8bit SI P
        public void OneRoundTrans_ShiftInvariant_53SPN(int[] table, int[] Ptable, int[] Stable5, int[] Stable3)
        {
            int xl, xr;//xl为5bit,xr为3bit
            int Sout = 0;
            for (int i = 0; i < len; i++)
            {
                //拆分
                xl = (table[i]) >> 3;
                xr = table[i] & 0x07;
                //S
                Sout = (Stable5[xl] << 3) | Stable3[xr];
                //P
                table[i] = Ptable[Sout];
            }
        }
        //5bit+3bit ShiftInvariant置换+8bit BitP
        public void OneRoundTrans_ShiftInvariant_53SBitP(int[] table, int[] Ptable, int[] Stable5, int[] Stable3)
        {
            int xl, xr;//xl为5bit,xr为3bit
            int Sout = 0;
            for (int i = 0; i < len; i++)
            {
                //拆分
                xl = (table[i]) >> 3;
                xr = table[i] & 0x07;
                //S
                Sout = (Stable5[xl] << 3) | Stable3[xr];
                //P
                table[i] = Ptable[Sout];
               
            }
        }
        //5bit+3bit ShiftInvariant置换+8bit BitP_固定常数C
        public void OneRoundTrans_ShiftInvariant_53SBitP_FixC(int[] table, int[] Ptable, int[] Stable5, int[] Stable3, int C)
        {
            int xl, xr;//xl为5bit,xr为3bit
            int Sout = 0;
            for (int i = 0; i < len; i++)
            {
                //拆分
                xl = (table[i]) >> 3;
                xr = table[i] & 0x07;
                //S
                Sout = (Stable5[xl] << 3) | Stable3[xr];
                //P
                table[i] = Ptable[Sout];
                table[i] = table[i] ^ C;
            }
        }
        //5-3卡方函数+拉线P
        public void OneRoundTrans_X2_BitP(int[] table, int[] Ptable)
        {
            int xl, xr;//xl为5bit,xr为3bit
            int xl1, xl2, xr1, xr2;//循环移位结果
            int y, yl, yr;
            int pout;
            for (int i = 0; i < len; i++)
            {
                //拆分
                xl = (table[i]) >> 3;
                xr = table[i] & 0x07;
                //卡方函数
                xl1 = ((xl << 1) | (xl >> 4))&0x1f;
                xl2 = ((xl << 2) | (xl >> 3))&0x1f;
                xr1 = ((xr << 1) | (xr >> 2))&0x07;
                xr2 = ((xr << 2) | (xr >> 1))&0x07;

                yl = xl ^ ((xl1 ^ 0x1f) & xl2);
                yr = xr ^ ((xr1 ^ 0x07) & xr2);
                //合并
                y = (yl << 3) | yr;
                //bitP
                pout = 0;
                for (int j = 0; j < 8; j++)
                    if (GetBit(y, Ptable[j]) == 1)
                        pout = SetBit(pout, j);
                //加常数
                table[i] = pout;
            }
        }
        //2分支，每个分支4比特
        public void OneRoundTrans_balanced(int[] table, long func,long func1)
        {
            int newX0, newX1;
            
            for (int i = 0; i < len; i++)
            {
                int X0 = table[i] >> 4;
                int X1 = table[i] & 0xf;
                int Fin = table[i] >> 4;
                int Fout = 0;
                for (int shift = 0; shift < 2; shift=shift+1)
                {
                    int mask = (0x01 << shift) - 1;
                    Fout = Fout | GetBit(func, ((Fin >> (shift)) | ((Fin & mask) << (4 - shift))));
                        Fout = Fout << 1;
                }
                for (int shift = 0; shift <2; shift=shift+2)
                {
                    int mask = (0x01 << shift) - 1;
                    Fout = Fout | GetBit(func1, ((Fin >> (shift)) | ((Fin & mask) << (4 - shift))));
                    if (shift != 1)
                        Fout = Fout << 1;
                }


                newX1 = X0 ;
                newX0 = X1^Fout;
                table[i] = ((newX0 << 4) | (newX1));
            }
        }

        //2分支，每个分支4比特
        public void OneRoundTrans_balanced1(int[] table, long func)
        {
            int newX0, newX1;

            for (int i = 0; i < len; i++)
            {
                int X0 = table[i] >> 4;
                int X1 = table[i] & 0xf;
                int Fin = table[i] >> 4;
                int Fout = 0;
                for (int shift = 0; shift < 4; shift++)
                {
                    int mask = (0x01 << shift) - 1;
                    Fout = Fout | GetBit(func, ((Fin >> (shift)) | ((Fin & mask) << (4 - shift))));
                    if (shift != 3)
                        Fout = Fout << 1;
                }

                newX1 = X0;
                newX0 = X1 ^ Fout;
                table[i] = ((newX0 << 4) | (newX1));
            }
        }
        //2分支，每个分支4比特
        public void OneRoundTrans_balanced_WithP(int[] table, long func,int[][] Pmatrix)
        {
            int newX0, newX1;

            for (int i = 0; i < len; i++)
            {
                int X0 = table[i] >> 4;
                int X1 = table[i] & 0xf;
                int Fin = table[i] >> 4;
                int[] Sout = new int[4];
                for (int shift = 0; shift < 4; shift++)
                {
                    int mask = (0x01 << shift) - 1;
                    Sout[shift]= GetBit(func, ((Fin >> (shift)) | ((Fin & mask) << (4 - shift))));
                }
                int Fout = 0;
                //矩阵乘
                for (int r = 0; r < 4; r++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        Fout = Fout ^ (Pmatrix[r][j] * Sout[j]);
                    }
                    if (r != 3)
                        Fout = Fout << 1;
                }


                newX1 = X0;
                newX0 = X1 ^ Fout;
                table[i] = ((newX0 << 4) | (newX1));
            }
        }

        //2分支，每个分支4比特
        public void OneRoundTrans_balanced_WithP_BitPerm(int[] table, long func,int[][] Pmatrix, int[] Ptable)
        {
            int newX0, newX1;

            for (int i = 0; i < len; i++)
            {
                int X0 = table[i] >> 4;
                int X1 = table[i] & 0xf;
                int Fin = table[i] >> 4;
                int[] Sout = new int[4];
                int temp = Fin;
                for (int shift = 0; shift < 4; shift++)
                {
                    Sout[shift] = GetBit(func, temp);
                    temp = Ptable[temp];
                }
                int Fout = 0;
                //矩阵乘
                for (int r = 0; r < 4; r++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        Fout = Fout ^ (Pmatrix[r][j] * Sout[j]);
                    }
                    if (r != 3)
                        Fout = Fout << 1;
                }


                newX1 = X0;
                newX0 = X1 ^ Fout;
                table[i] = ((newX0 << 4) | (newX1));
            }
        }
        //
        public void OneRoundTrans_Type2(int[] table, long[] func)
        {
            int newX0, newX1, newX2, newX3;

            for (int i = 0; i < len; i++)
            {
               
                int x0 = (GetBit(table[i], 1) << 1) | GetBit(table[i], 0);
                int x1 = (GetBit(table[i], 3) << 1) | GetBit(table[i], 2);
                int x2 = (GetBit(table[i], 5) << 1) | GetBit(table[i], 4);
                int x3 = (GetBit(table[i], 7) << 1) | GetBit(table[i], 6);
                int F1out = (GetBit(func[3], x3) << 1) | GetBit(func[2], x3);
                int F0out = (GetBit(func[1], x1) << 1) | GetBit(func[0], x1);

                newX3 = F1out ^ x2;
                newX2 = x1;
                newX1 = F0out ^ x0;
                newX0 = x3;
                table[i] = ((newX3 << 6) | (newX2 << 4) | (newX1 << 2) | newX0);
            }
        }

        public void OneRoundTrans_Type2_AddConstant(int[] table, long[] func,int rc)
        {
            int newX0, newX1, newX2, newX3;

            for (int i = 0; i < len; i++)
            {

                int x0 = (GetBit(table[i], 1) << 1) | GetBit(table[i], 0);
                int x1 = (GetBit(table[i], 3) << 1) | GetBit(table[i], 2);
                int x2 = (GetBit(table[i], 5) << 1) | GetBit(table[i], 4);
                int x3 = (GetBit(table[i], 7) << 1) | GetBit(table[i], 6);
                int F1out = (GetBit(func[3], x3 ^ ((rc >> 2) & 0x3)) << 1) | GetBit(func[2], x3 ^ ((rc >> 2) & 0x3));
                int F0out = (GetBit(func[1], x1 ^ (rc & 0x3)) << 1) | GetBit(func[0], x1 ^ (rc & 0x3));

                newX3 = F1out ^ x2;
                newX2 = x1 ;
                newX1 = F0out ^ x0 ;
                newX0 = x3 ;
                table[i] = ((newX3 << 6) | (newX2 << 4) | (newX1 << 2) | newX0);
            }
        }
        //
        public void OneRoundTrans_Type3(int[] table, long[] func)
        {
            int newX0, newX1, newX2, newX3;

            for (int i = 0; i < len; i++)
            {

                int x0 = (GetBit(table[i], 1) << 1) | GetBit(table[i], 0);
                int x1 = (GetBit(table[i], 3) << 1) | GetBit(table[i], 2);
                int x2 = (GetBit(table[i], 5) << 1) | GetBit(table[i], 4);
                int x3 = (GetBit(table[i], 7) << 1) | GetBit(table[i], 6);
                int F2out = (GetBit(func[5], x3) << 1) | GetBit(func[4], x3);
                int F1out = (GetBit(func[3], x2) << 1) | GetBit(func[2], x2);
                int F0out = (GetBit(func[1], x1) << 1) | GetBit(func[0], x1);
                newX3 = F2out ^ x2;
                newX2 = F1out ^ x1;
                newX1 = F0out ^ x0;
                newX0 = x3;
                table[i] = ((newX3 << 6) | (newX2 << 4) | (newX1 << 2) | newX0);
            }
        }

        public void OneRoundTrans_Type3_AddConstant(int[] table, long[] func,int rc)
        {
            int newX0, newX1, newX2, newX3;

            for (int i = 0; i < len; i++)
            {

                int x0 = (GetBit(table[i], 1) << 1) | GetBit(table[i], 0);
                int x1 = (GetBit(table[i], 3) << 1) | GetBit(table[i], 2);
                int x2 = (GetBit(table[i], 5) << 1) | GetBit(table[i], 4);
                int x3 = (GetBit(table[i], 7) << 1) | GetBit(table[i], 6);
                int F2out = (GetBit(func[5], x3 ^ ((rc >> 4) & 0x3)) << 1) | GetBit(func[4], x3 ^ ((rc >> 4) & 0x3));
                int F1out = (GetBit(func[3], x2 ^ ((rc >> 2) & 0x3)) << 1) | GetBit(func[2], x2 ^ ((rc >> 2) & 0x3));
                int F0out = (GetBit(func[1], x1 ^ ((rc >> 0) & 0x3)) << 1) | GetBit(func[0], x1 ^ ((rc >> 0) & 0x3));
                newX3 = F2out ^ x2;
                newX2 = F1out ^ x1;
                newX1 = F0out ^ x0;
                newX0 = x3;
                table[i] = ((newX3 << 6) | (newX2 << 4) | (newX1 << 2) | newX0);
            }
        }
        public void OneRoundTrans_balanced_BitPerm(int[] table, long func,int[] Ptable)
        {
            int newX0, newX1;

            for (int i = 0; i < len; i++)
            {
                int X0 = table[i] >> 4;
                int X1 = table[i] & 0xf;
                int Fin = table[i] >> 4;
                int Fout = 0;
                int temp=Fin;
                for (int j = 0;j < 4; j++)
                {
                    Fout = Fout | GetBit(func, temp);
                    temp=Ptable[temp];
                    if (j != 3)
                        Fout = Fout << 1;
                }

                newX1 = X0;
                newX0 = X1 ^ Fout;
                table[i] = ((newX0 << 4) | (newX1));
            }
        }
        
        //2分支，每个分支4比特
        public void OneRoundTrans_balanced_PS(int[] table, long func, int[][] Pmatrix)
        {
            int newX0, newX1;

            for (int i = 0; i < len; i++)
            {
                int X0 = table[i] >> 4;
                int X1 = table[i] & 0xf;
                int Fin = table[i] >> 4;
                int[] Finb = new int[4];
                for (int j = 0; j < 4; j++)
                    Finb[j] = GetBit(Fin, j);
                //P
                int Pout =0;
                //矩阵乘
                for (int r = 0; r < 4; r++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                      Pout  = Pout ^ (Pmatrix[r][j] * Finb[j]);
                    }
                    if (r != 3)
                        Pout = Pout << 1;
                }
                int Fout = 0;

                for (int shift = 0; shift < 4; shift++)
                {
                    int mask = (0x01 << shift) - 1;
                    Fout = GetBit(func, ((Pout >> (shift)) | ((Pout & mask) << (4 - shift))));
                    if (shift != 3)
                        Fout = Fout << 1;
                }
                

                newX1 = X0;
                newX0 = X1 ^ Fout;
                table[i] = ((newX0 << 4) | (newX1));
            }
        }
        /// <summary>
        /// 寻找最优S盒：
        /// table表示S盒的真值表 (6进2出4分支)
        /// </summary>
        /// <param name="table"></param>
        /// <param name="size"></param>
        /// <param name="round"></param>
        /// <param name="optimalDiff"></param>
        /// <param name="optimalNonlinear"></param>
        /// <returns></returns>
        public int[] SearchOptimal(string filename, string scriptfilename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];

            //fsScript = new FileStream(scriptfilename, FileMode.Create);
            //swScript = new StreamWriter(fsScript);

            //写Script头
            //swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            //swScript.WriteLine("set synthetic_library generic.sdb");
            //swScript.WriteLine("set target_library \"fast.db slow.db\"");
            //swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");



            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);

            //每轮的布尔函数都是一样的
            int varNum =4;
            int[] ANF = new int[Power(2, varNum)];
            int terms = 1 + varNum + varNum * (varNum - 1) / 2;//ANF项数量，1+6+6C2
            int limit =0x1<<terms;
            for (int num = 0; num < limit; num++)
            {
               // if (num % 65536 == 0)
                //    System.Console.WriteLine("num={0:x}", num);
                //从num中获取ANF
                for (int i = 0; i < terms; i++)
                {
                    if (GetBit(num, i) == 1)
                        ANF[RankTable[i]] = 1;
                    else
                        ANF[RankTable[i]] = 0;
                }
                //从ANF获得真值表
                long func = InvMoebiusTrans(varNum, ANF);

                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                //根据num的值，确定轮函数
                for (int i = 0; i < round; i++)
                {
                    item[i] = func;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans(table, item[i]);
                }

                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);
                
                if((num&0xffff)==0)
                System.Console.WriteLine("diff={0},num={1:x}", diff, num);
                
                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数真值表={0:x}",func);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();
                        //swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/F_I6O1_D2_" + Convert.ToString(func,16) + ".v\"}");
                        //swScript.WriteLine("compile -exact_map");
                        //swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/areareport_" + num + ".txt\"");
                        //swScript.WriteLine("remove_design -designs");
                        //return table;
                    }
                    
                }
            }
            System.Console.WriteLine("共有{0}个可行解", count);
            sw.Close();
            fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            return errorTable;
            //释放item
        }
        //检验是否为置换
        public bool CheckPermutaion(int[] table)
        {
            int[] temp = new int[table.Length];
            Array.Copy(table, temp, table.Length);
            Array.Sort(temp);
            for (int i = 0; i < table.Length; i++)
                if (temp[i] != i)
                    return false;
            return true;
        }
        //小数量的阶乘函数
        public int factorial(int x)
        {
            int result = 1;
            for (int i = 1; i <= x; i++)
                result *= i;
            return result;
        }
        //从序号换算8比特置换表
        //
        public int[] GetBitPerm8(int no)
        {
            int[] BitP = new int[8];
            int[] rank = new int[8];
            bool[] Used = new bool[8];
            //换算rank
            for (int i = 0; i < 8; i++)
            {
                rank[i] = no / factorial(7 - i);
                no = no % factorial(7 - i);
                //根据rank[i]找相应的比特
                int t = 0;
                int count = 0;
                while (count < rank[i] || (Used[t]))//停止条件：等于rank且没有用过
                {
                    if (Used[t])//用过的直接跳过
                    {
                        t++;
                        continue;
                    }
                    //没用过的
                    if (count < rank[i])
                        count++;
                    //若已经到了rank但是当前t用过，则直接跳过
                    t++;

                }
                if (t >= 8)
                    System.Console.WriteLine("Error!");
                Used[t] = true;
                BitP[i] = t;
            }
            return BitP;
        }
        /// <summary>
        /// 寻找最优S盒：SP结构，5-3 X2函数做S，P为比特置换，加固定常数
        /// 复杂度：8！*2^8
        /// </summary>
        /// <param name="table"></param>
        /// <param name="size"></param>
        /// <param name="round"></param>
        /// <param name="optimalDiff"></param>
        /// <param name="optimalNonlinear"></param>
        /// <returns></returns>
        public int[] SearchOptimal_X2_BitP_FixC(string filename, string scriptfilename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //fsScript = new FileStream(scriptfilename, FileMode.Create);
            //swScript = new StreamWriter(fsScript);

            //写Script头
            //swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            //swScript.WriteLine("set synthetic_library generic.sdb");
            //swScript.WriteLine("set target_library \"fast.db slow.db\"");
            //swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");



            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);

            for (int num = 0; num < 40320*256; num++)
            {
                if (num % 65536 == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, num/(40320 * 2.56));
                }
                //从num中获取常数C和比特置换P
                int P = num >> 8;
                int C = num & 0xff;
                //从P中获取置换表
                int[] Ptable = GetBitPerm8(P);

                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_X2_BitP_FixC(table, Ptable,C);
                }
                //if (!CheckPermutaion(table))
                //{
                //    System.Console.WriteLine("None permuataion!");
                //    return errorTable;
                //};
                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    mindiff = diff;

                if ((num & 0xffff) == 0)
                    System.Console.WriteLine("diff={0},num={1:x}", diff, num);

                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    //System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear > maxnL)
                        maxnL = nonLinear;
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数={0:x}", num);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();
                        //swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/F_I6O1_D2_" + Convert.ToString(func,16) + ".v\"}");
                        //swScript.WriteLine("compile -exact_map");
                        //swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/areareport_" + num + ".txt\"");
                        //swScript.WriteLine("remove_design -designs");
                        //return table;
                    }

                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            return errorTable;
        }
        public int[] SearchOptimal_X2_BitP(string filename, string scriptfilename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //fsScript = new FileStream(scriptfilename, FileMode.Create);
            //swScript = new StreamWriter(fsScript);

            //写Script头
            //swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            //swScript.WriteLine("set synthetic_library generic.sdb");
            //swScript.WriteLine("set target_library \"fast.db slow.db\"");
            //swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");



            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);

            for (int num = 0; num < 40320; num++)
            {
                if (num % 65536 == 0)
                    System.Console.WriteLine("num={0:x}", num);
                //从num中获取常数C和比特置换P
                int P = num ;

                //从P中获取置换表
                int[] Ptable = GetBitPerm8(P);

                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_X2_BitP(table, Ptable);
                    //if (!CheckPermutaion(table))
                    //{
                    //    System.Console.WriteLine("None permuataion!");
                    //    return errorTable;
                    //};
                }
                if (!CheckPermutaion(table))
                {
                    System.Console.WriteLine("None permuataion!");
                    return errorTable;
                };
                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    mindiff = diff;

                if ((num & 0xffff) == 0)
                    System.Console.WriteLine("diff={0},num={1:x}", diff, num);

                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear > maxnL)
                        maxnL = nonLinear;
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数={0:x}", num);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();
                        //swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/F_I6O1_D2_" + Convert.ToString(func,16) + ".v\"}");
                        //swScript.WriteLine("compile -exact_map");
                        //swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/areareport_" + num + ".txt\"");
                        //swScript.WriteLine("remove_design -designs");
                        //return table;
                    }

                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            return errorTable;
            //释放item
        }
        public int[] SearchOptimal_Type2(string filename, string scriptfilename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff=256;
            int maxnL = 0;

            //fsScript = new FileStream(scriptfilename, FileMode.Create);
            //swScript = new StreamWriter(fsScript);

            //写Script头
            //swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            //swScript.WriteLine("set synthetic_library generic.sdb");
            //swScript.WriteLine("set target_library \"fast.db slow.db\"");
            //swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");



            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);

            //每轮的布尔函数都是一样的
          
            int limit = 0x1 << 16;
             long[] func = new long[4];
            for (int num = 0; num < limit; num++)
            {
                
                //获得真值表
                func[0] = num & 0xf;
                func[1] = (num >> 4) & 0xf;
                func[2] = (num >> 8) & 0xf;
                func[3] = (num >> 12) & 0xf;

                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                //根据num的值，确定轮函数

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_Type2(table, func);
                }

                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    mindiff = diff;

                if ((num & 0xff) == 0)
                    System.Console.WriteLine("diff={0},num={1:x}", diff, num);

                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear > maxnL)
                        maxnL = nonLinear;
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数真值表={0:x}", func);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();
                        //swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/F_I6O1_D2_" + Convert.ToString(func,16) + ".v\"}");
                        //swScript.WriteLine("compile -exact_map");
                        //swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/areareport_" + num + ".txt\"");
                        //swScript.WriteLine("remove_design -designs");
                        //return table;
                    }

                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            return errorTable;
            //释放item
        }

        public int[] SearchOptimal_Type2_AddConstant(string filename, string scriptfilename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            int Constant = 1;
            int rc=0;
            //fsScript = new FileStream(scriptfilename, FileMode.Create);
            //swScript = new StreamWriter(fsScript);

            //写Script头
            //swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            //swScript.WriteLine("set synthetic_library generic.sdb");
            //swScript.WriteLine("set target_library \"fast.db slow.db\"");
            //swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");



            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);

            //每轮的布尔函数都是一样的

            int limit = 0x1 << 16;
            long[] func = new long[4];
            for (int num = 0; num < limit; num++)
            {

                //获得真值表
                func[0] = num & 0xf;
                func[1] = (num >> 4) & 0xf;
                func[2] = (num >> 8) & 0xf;
                func[3] = (num >> 12) & 0xf;

                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                //根据num的值，确定轮函数

                rc = Constant;

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_Type2_AddConstant(table, func,rc);
                    rc = (rc >> 1) | (rc & 0x1) << 3;
                }

                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    mindiff = diff;

                if ((num & 0xff) == 0)
                    System.Console.WriteLine("diff={0},num={1:x}", diff, num);

                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear > maxnL)
                        maxnL = nonLinear;
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数真值表={0:x}", func);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();
                        //swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/F_I6O1_D2_" + Convert.ToString(func,16) + ".v\"}");
                        //swScript.WriteLine("compile -exact_map");
                        //swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/areareport_" + num + ".txt\"");
                        //swScript.WriteLine("remove_design -designs");
                        //return table;
                    }

                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            return errorTable;
            //释放item
        }

        public int[] SearchOptimal_Type3(string filename, string scriptfilename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;

            //fsScript = new FileStream(scriptfilename, FileMode.Create);
            //swScript = new StreamWriter(fsScript);

            //写Script头
            //swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            //swScript.WriteLine("set synthetic_library generic.sdb");
            //swScript.WriteLine("set target_library \"fast.db slow.db\"");
            //swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");



            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);

            //每轮的布尔函数都是一样的

            int limit = 0x1 << 24;
            long[] func = new long[6];
            for (int num = 0; num < limit; num++)
            {
                if ((num & 0xffff) == 0)
                    System.Console.WriteLine("mindiff={0},num={1:x}", mindiff, num);
                //获得真值表
                func[0] = num & 0xf;
                func[1] = (num >> 4) & 0xf;
                func[2] = (num >> 8) & 0xf;
                func[3] = (num >> 12) & 0xf;
                func[4] = (num >> 16) & 0xf;
                func[5] = (num >> 20) & 0xf;
                bool remove=false;
                for (int i = 0; i < 6; i++)
                {
                    if (func[i] == 0 || func[i] == 15)
                    {
                        remove = true;
                        break;
                    }
                }
                if (remove)
                    continue;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                //根据num的值，确定轮函数

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_Type3(table, func);
                }

                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    mindiff = diff;



                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear > maxnL)
                        maxnL = nonLinear;
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数真值表={0:x}", func);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();
                        //swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/F_I6O1_D2_" + Convert.ToString(func,16) + ".v\"}");
                        //swScript.WriteLine("compile -exact_map");
                        //swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/areareport_" + num + ".txt\"");
                        //swScript.WriteLine("remove_design -designs");
                        //return table;
                    }

                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            return errorTable;
            //释放item
        }

        public int[] SearchOptimal_Type3_AddConstant(string filename, string scriptfilename, int optimalDiff, int optimalNonlinear,int Constant)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int rc = 0;
            int mindiff = 256;
            int maxnL = 0;

            //fsScript = new FileStream(scriptfilename, FileMode.Create);
            //swScript = new StreamWriter(fsScript);

            //写Script头
            //swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            //swScript.WriteLine("set synthetic_library generic.sdb");
            //swScript.WriteLine("set target_library \"fast.db slow.db\"");
            //swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");



            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);

            //每轮的布尔函数都是一样的

            int limit = 0x1 << 24;
            long[] func = new long[6];
            for (int num = 0; num < limit; num++)
            {
                if ((num & 0xffff) == 0)
                    System.Console.WriteLine("mindiff={0},num={1:x}", mindiff, num);
                //获得真值表
                func[0] = num & 0xf;
                func[1] = (num >> 4) & 0xf;
                func[2] = (num >> 8) & 0xf;
                func[3] = (num >> 12) & 0xf;
                func[4] = (num >> 16) & 0xf;
                func[5] = (num >> 20) & 0xf;
                bool remove = false;
                for (int i = 0; i < 6; i++)
                {
                    if (func[i] == 0 || func[i] == 15)
                    {
                        remove = true;
                        break;
                    }
                }
                if (remove)
                    continue;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }
                rc = Constant;
                //根据num的值，确定轮函数

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_Type3_AddConstant(table, func,rc);
                    rc = (rc >> 1) | ((rc & 0x1) << 5);
                }

                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    mindiff = diff;



                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear > maxnL)
                        maxnL = nonLinear;
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数真值表={0:x}", func);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();
                        //swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/F_I6O1_D2_" + Convert.ToString(func,16) + ".v\"}");
                        //swScript.WriteLine("compile -exact_map");
                        //swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/areareport_" + num + ".txt\"");
                        //swScript.WriteLine("remove_design -designs");
                        //return table;
                    }

                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            return errorTable;
            //释放item
        }
        //二元阵求逆
        #region
        //
        //函数GJ_elimination，对0-1矩阵进行高斯消元
        public void GJ_elimination(int[][] M, int rownum, int colnum)
        {
            int i, j, k, flag, temp;

            int rank = 0;

            int startcol, currrow, postrow=0;

            startcol = 0;   //从矩阵的第[0,0]位置开始实施消元算法

            for (i = 0; i < rownum; i++)   //i代表第i次消元过程
            {
                currrow = i;

                flag = 0;

                //寻找第i次高斯消元的主元素所在的行列位置
                for (j = startcol; j < colnum; j++)
                {
                    for (k = currrow; k < rownum; k++)
                    {
                        if (M[k][j] != 0)
                        {
                            startcol = j + 1;

                            postrow = k;

                            flag = 1;

                            (rank)++;

                            break;
                        }
                    }

                    if (flag == 1)
                    {
                        break;
                    }
                }

                if (flag == 0)
                    break;   //表示高斯消元其实已经完成，下面的行中元素都为0；跳出i循环；

                else
                {
                    if (currrow != postrow)   //如果currrow!=postrow,则交换两行的元素；
                    {
                        for (j = startcol - 1; j < colnum; j++)
                        {
                            temp = M[currrow][j];
                            M[currrow][j] = M[postrow][j];
                            M[postrow][j] = temp;
                        }
                    }

                    for (j = 0; j < rownum; j++)   //用currrow对其它的行进行高斯消元
                    {
                        if (j != currrow && M[j][startcol - 1] != 0)
                        {
                            for (k = startcol - 1; k < colnum; k++)
                                M[j][k] = (M[j][k] ^ M[currrow][k]);
                        }
                    }
                }
            }
        }

        //函数ReverseM，求“可逆矩阵”的逆
        public int ReverseM(int[][] M, int dim) //dim == rownum
        {
            int i, j;

            //初始化一个新矩阵等于[M,I]
            int[][] tempM = new int[dim][];

            for (i = 0; i < dim; i++)
            {
                tempM[i] = new int[dim*2];
            }

            for (i = 0; i < dim; i++)
            {
                for (j = 0; j < dim; j++)
                {
                    tempM[i][j] = M[i][j];
                }
            }

            for (i = 0; i < dim; i++)
            {
                for (j = dim; j < 2 * dim; j++)
                {
                    if (j == (dim + i))
                        tempM[i][j] = 1;

                    else
                        tempM[i][j] = 0;
                }
            }

            GJ_elimination(tempM, dim, 2 * dim);

            if (tempM[dim - 1][dim - 1] == 0)
                return 0;

            else
            {
              
                return 1;
            }
        }
        #endregion

        public int[] CreatePermTable(int l, int[] Perm)
        {
            int len = (int)Math.Pow(2, l);
            int[] table=new int[len];
            for (int i = 0; i < len; i++)
            {
                table[i] = 0;
                for (int j = 0; j < l; j++)
                    if (GetBit(i, j) != 0)
                       table[i]= SetBit(table[i], Perm[j]);
            }

            //检查
            //Array.Sort(table);
            //for (int i = 0; i < len - 1; i++)
            //    if (table[i] == table[i + 1])
            //        return null;
            //检查结束
            return table;
        }
        public int[] SearchOptimal_WithP_BitPerm(string filename, int optimalDiff, int optimalNonlinear)
        {
            //重新生成ranktable
            int mindiff = 256;
            int maxNl = 0;
            int nsize = size / 4 * 2;
            int m = 1 + nsize + nsize * (nsize - 1) / 2 ;
            RankTable = new int[m];
            RankTable[0] = 0;
            HammingWeight hw = new HammingWeight(nsize);
            int i = 1;
            while (i < m)
            {
                //  hw.PrintState();
                RankTable[i] = hw.ReturnNum();
                hw.HwNext();
                i++;
            }


            //-----------


            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1}", round, size);

            //每轮的布尔函数都是一样的
            int varNum = 4;
            int[] ANF = new int[Power(2, varNum)];
            int terms = m+4*4;//ANF项数量，1+4+4C2
            int limit = 0x1 << terms;
            long num=0;
            for (int pind =0; pind < (0x1 << 8); pind++)
            {

                //生成比特置换组
                int[] Perm = new int[varNum];
                int[] Perm_Copy = new int[varNum];
                int temp = pind;
                for (i = 0; i < varNum; i++)
                {
                    Perm[i] = (int)(temp & 0x03);
                    Perm_Copy[i] = Perm[i];
                    temp = temp >> 2;
                }
                Array.Sort(Perm_Copy);
                //检查是否有重复
                bool Pvalid = true;
                for (i = 0; i < varNum - 1; i++)
                    if (Perm_Copy[i] == Perm_Copy[i + 1])
                    {
                        Pvalid = false;
                        break;
                    }
                if (!Pvalid)
                    continue;
                //生成P对应的置换表
                int[] Ptable = CreatePermTable(varNum, Perm);

                System.Console.WriteLine("mindiff={0},pind={1:x}", mindiff,pind);
                for (long pnum = 0; pnum < (0x1 << 16); pnum++)
                {
                    // if (num % 65536 == 0)
                    //    System.Console.WriteLine("num={0:x}", num);

                    int[][] Pmatrix = new int[varNum][];
                    for (i = 0; i < varNum; i++)
                    {
                        Pmatrix[i] = new int[varNum];
                        for (int j = 0; j < varNum; j++)
                            Pmatrix[i][j] = GetBit(pnum, i * varNum + j);
                    }
                    //矩阵不可逆时跳过
                    if (ReverseM(Pmatrix, varNum) == 0)
                        continue;
                    System.Console.WriteLine("mindiff={0},pnum={1:x}", mindiff, pnum);
                    for (long snum = 0; snum < (0x1 << m); snum++)
                    {
                        num = (snum << 24) |(pnum<<8)|((long)pind);
                        //从num中获取ANF
                        for (i = 0; i < m; i++)
                        {
                            if (GetBit(snum, i) == 1)
                                ANF[RankTable[i]] = 1;
                            else
                                ANF[RankTable[i]] = 0;
                        }



                        //从ANF获得真值表
                        long func = InvMoebiusTrans(varNum, ANF);

                        //table初始化为恒等变换
                        for (i = 0; i < len; i++)
                        {
                            table[i] = i;
                        }

                        //根据num的值，确定轮函数
                        for (i = 0; i < round; i++)
                        {
                            item[i] = func;
                        }

                        for (i = 0; i < round; i++)
                        {
                            OneRoundTrans_balanced_WithP_BitPerm(table, item[i], Pmatrix, Ptable);
                        }

                        //求解最终S盒的差分均匀度、非线性度
                        int diff = DiffUniform(table, size);
                        if (diff < mindiff)
                            mindiff = diff;
                        //if ((num & 0xffff000000) == 0)
                       //     System.Console.WriteLine("diff={0},num={1:x}", diff, num);

                        if (diff <= optimalDiff)
                        {
                            int nonLinear = Nonlinear(table, size);
                            if (nonLinear > maxNl)
                                maxNl = nonLinear;
                            System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                            if (nonLinear >= optimalNonlinear)
                            {
                                System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                                count++;
                                sw.WriteLine("******************");
                                sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                                sw.WriteLine("输出表：");
                                for (int j = 0; j < 256; j++)
                                    sw.Write("{0:x2}\t", table[j]);
                                sw.WriteLine("\n");
                                sw.WriteLine("内部函数真值表={0:x}", func);
                                sw.WriteLine("******************");
                                sw.WriteLine("\n");
                                sw.WriteLine("\n");
                                sw.WriteLine("\n");
                                sw.Flush();

                            }

                        }
                    }
                }
            }
            System.Console.WriteLine("共有{0}个可行解", count);
            System.Console.WriteLine("最小差分={0},最大非线性={1}", mindiff,maxNl);
            sw.WriteLine("最小差分={0},最大非线性={1}", mindiff, maxNl);
            sw.Close();
            fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            return errorTable;
            //释放item
        }

        public int[] SearchOptimal_BitPerm(string filename, int optimalDiff, int optimalNonlinear)
        {
            //重新生成ranktable
            int mindiff = 256;
            int maxNl = 0;
            int nsize = size / 4 * 2;
            int m = 1 + nsize + nsize * (nsize - 1) / 2;
            RankTable = new int[m];
            RankTable[0] = 0;
            HammingWeight hw = new HammingWeight(nsize);
            int i = 1;
            while (i < m)
            {
                //  hw.PrintState();
                RankTable[i] = hw.ReturnNum();
                hw.HwNext();
                i++;
            }


            //-----------


            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1}", round, size);

            //每轮的布尔函数都是一样的
            int varNum = 4;
            int[] ANF = new int[Power(2, varNum)];
            int terms = m + 8;//ANF项数量，1+4+4C2
            int limit = 0x1 << terms;
            long num = 0;
            for (long pnum = 0; pnum < (0x1 <<8); pnum++)
            {
                // if (num % 65536 == 0)
                //    System.Console.WriteLine("num={0:x}", num);
                //生成比特置换组
                int[] Perm = new int[varNum];
                int[] Perm_Copy = new int[varNum];
                long temp=pnum;
                for (i = 0; i < varNum; i++)
                {
                    Perm[i] = (int)(temp&0x03);
                    Perm_Copy[i] = Perm[i];
                    temp = temp >> 2;
                }
                Array.Sort(Perm_Copy);
                //检查是否有重复
                bool Pvalid=true;
                for(i=0;i<varNum-1;i++)
                    if (Perm_Copy[i] == Perm_Copy[i + 1])
                    {
                        Pvalid = false;
                        break;
                    }
                if (!Pvalid)
                    continue;
                //生成P对应的置换表
                int[] Ptable = CreatePermTable(varNum, Perm);

                for (long snum = 0; snum < (0x1 << m); snum++)
                {
                    num = (snum << 8) | pnum;
                    //从num中获取ANF
                    for (i = 0; i < m; i++)
                    {
                        if (GetBit(snum, i) == 1)
                            ANF[RankTable[i]] = 1;
                        else
                            ANF[RankTable[i]] = 0;
                    }



                    //从ANF获得真值表
                    long func = InvMoebiusTrans(varNum, ANF);

                    //table初始化为恒等变换
                    for (i = 0; i < len; i++)
                    {
                        table[i] = i;
                    }

                    //根据num的值，确定轮函数
                    for (i = 0; i < round; i++)
                    {
                        item[i] = func;
                    }

                    for (i = 0; i < round; i++)
                    {
                        OneRoundTrans_balanced_BitPerm(table, item[i],Ptable);
                        //OneRoundTrans_balanced1(table, item[i]);
                    }

                    //求解最终S盒的差分均匀度、非线性度
                    int diff = DiffUniform(table, size);
                    if (diff < mindiff)
                        mindiff = diff;
                    if ((snum & 0xffff) == 0)
                        System.Console.WriteLine("diff={0},num={1:x}", diff, num);

                    if (diff <= optimalDiff)
                    {
                        int nonLinear = Nonlinear(table, size);
                        if (nonLinear > maxNl)
                            maxNl = nonLinear;
                        System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        if (nonLinear >= optimalNonlinear)
                        {
                            System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                            count++;
                            sw.WriteLine("******************");
                            sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                            sw.WriteLine("输出表：");
                            for (int j = 0; j < 256; j++)
                                sw.Write("{0:x2}\t", table[j]);
                            sw.WriteLine("\n");
                            sw.WriteLine("内部函数真值表={0:x}", func);
                            sw.WriteLine("******************");
                            sw.WriteLine("\n");
                            sw.WriteLine("\n");
                            sw.WriteLine("\n");
                            sw.Flush();

                        }

                    }
                }
            }
            System.Console.WriteLine("共有{0}个可行解", count);
            System.Console.WriteLine("最小差分={0},最大非线性={1}", mindiff, maxNl);
            sw.WriteLine("最小差分={0},最大非线性={1}", mindiff, maxNl);
            sw.Close();
            fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            return errorTable;
            //释放item
        }
        //计算HW
        public int HW(int x, int len)
        {
            int c = 0;
            for(int i=0;i<len;i++)
                if(GetBit(x,i)!=0)
                   c++;
            return c;
        }
        //根据ANF算最大次数,v表示变量个数
        public int MaxDegree(int[] ANF,int v)
        {
            int maxd = 0;
            for (int i = 0; i < ANF.Length; i++)
            {
                if (ANF[i] == 1)
                {
                    int d = HW(i, v);
                    if (d > maxd)
                        maxd = d;
                }
            }
            return maxd;
        }
        public int[] SearchOptimal_balanced1(string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            int[] item = new int[round];
             int mindiff = 256;
            int maxNl = 0;
            int count1 = 0;

            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);


            //每轮的布尔函数都是一样的
            int varNum = 4;
            int[] ANF = new int[Power(2, varNum)];

            for (int num = 0; num < 256*256; num++)
            {
                MoebiusTrans(varNum, num, ANF);//真值表换ANF
                int maxd = MaxDegree(ANF, varNum);
                if (maxd > 2)
                    continue;//次数大于2的不算
                count1++;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                //根据num的值，确定轮函数
                for (int i = 0; i < round; i++)
                {
                    item[i] = num;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_balanced1(table, item[i]);
                }

                 //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);
                if (diff < mindiff)
                    mindiff = diff;
                if ((num & 0xffff) == 0)
                    System.Console.WriteLine("diff={0},num={1:x}", diff, num);

                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    if (nonLinear > maxNl)
                        maxNl = nonLinear;
                    System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数真值表={0:x}", num);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();

                    }

                }
            }
            System.Console.WriteLine("共有{0}个可行解", count);
            System.Console.WriteLine("最小差分={0},最大非线性={1}", mindiff,maxNl);
            
            System.Console.WriteLine("共有{0}个次数<=2的ANF", count1);
            sw.Close();
            fs.Close();



            return errorTable;
            //释放item
        }

        //从DC的areareport中提取areacost
        public double ReadAreaCost(string filename)
        {
            FileStream fs = new FileStream(filename, FileMode.Open);
            StreamReader sr = new StreamReader(fs);
            string line = sr.ReadLine();
            double cost = 0;
            while (line != null)
            {
                if (line.StartsWith("Total cell area:"))
                {
                    //int e=line.LastIndexOf(' ');
                    int s = line.Substring(0, line.Length - 1).LastIndexOf(' ') + 1;
                    cost = Double.Parse(line.Substring(s));
                };
                line = sr.ReadLine();
            }
            sr.Close();
            fs.Close();
            return cost;
        }
        //最后一个"_"和"."之间的部分为num
        public void GetCost_GE(string path, string CostFile)
        {
            FileStream fs = new FileStream(CostFile, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            DirectoryInfo TheFolder = new DirectoryInfo(path);
            //遍历文件
            foreach (FileInfo NextFile in TheFolder.GetFiles())
            {
                double cost = ReadAreaCost(NextFile.FullName) / 9.9792;
                string fname = NextFile.Name;
                int s = fname.LastIndexOf('_') + 1;
                int e = fname.LastIndexOf('.');
                long num = Convert.ToInt64(fname.Substring(s, e - s),16);
                sw.WriteLine("{0:x}\t{1}", num, cost);
            }
            sw.Close();
            fs.Close();
        }

        //根据SI规则构造完整真值表
        public int[] ShiftInvariantTT(byte[] TT, int bitlen)
        {
            int len = 0x1 << bitlen;
            int[] table=new int[len];
            for (int x = 0; x < len; x++)
            {
                table[x]=0;
                int xt = x;
                for (int bit = 0; bit < bitlen; bit++)
                {
                    if(GetBit(TT, xt)==1)
                        table[x] = SetBit(table[x], bit); 
                    
                     xt = (xt >> 1) | ((xt & 0x1) << (bitlen - 1));
                }
            }
            return table;
        }

        //两个4bit SI Sbox+8bit SI 线性变换
        public int[][] ReadShiftInvariantTT(string TTfile, int bitlen)
        {
            int tablelen = 0x1 << bitlen;
            int bytelen = tablelen / 8;
            int num = 0;
            FileStream fs = new FileStream(TTfile, FileMode.Open);
            BinaryReader br = new BinaryReader(fs);
            //先计算个数
            while (br.BaseStream.Position != br.BaseStream.Length)
            {
                br.ReadBytes(bytelen);
                num++;
            }
            //根据个数申请空间
            int[][] table=new int[num][];
            //读写指针还原
            br.BaseStream.Seek(0, SeekOrigin.Begin);
            for (int i = 0; i < num; i++)
            {
                table[i] = ShiftInvariantTT(br.ReadBytes(bytelen),bitlen);
            }
            br.Close();
            fs.Close();
            return table;
        }
        //两个4bit SI Sbox+8bit SI 线性变换
        public int[] SearchOptimal_ShiftInvariant_4SPN(string Sbinfile, string Pbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的4bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT(Sbinfile, 4);
            //读取所有有直接TI的8bit SI 线性置换
            int[][] Ptable = ReadShiftInvariantTT(Pbinfile, 8);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            long length = Stable.Length * Stable.Length * Ptable.Length;
            for (long num = Stable.Length * Stable.Length; num < length; num++)
            {
                if ((num&0xffff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0*num / length);
                }
                //从num中获取P的序号
                int Pno = (int)(num / (Stable.Length * Stable.Length));
                int Sno = (int)(num % (Stable.Length * Stable.Length));
                int Sno1 = Sno / Stable.Length;
                int Sno2 = Sno % Stable.Length;
                //从P中获取置换表
                //int[] Ptable = GetBitPerm8(P);

                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_ShiftInvariant_SPN(table, Ptable[Pno], Stable[Sno1],Stable[Sno2]);
                }
                //if (!CheckPermutaion(table))
                //{
                //    System.Console.WriteLine("None permuataion!");
                //    return errorTable;
                //};
                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    mindiff = diff;

                if ((num & 0xffff) == 0)
                    System.Console.WriteLine("diff={0},num={1:x}", diff, num);

                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    //System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear > maxnL)
                        maxnL = nonLinear;
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数={0:x}", num);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();

                    }

                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();


            return errorTable;
        }
        //5bit/3bit SI Sbox+8bit SI 线性变换
        public int[] SearchOptimal_ShiftInvariant_53SPN(string S5binfile, string S3binfile, string Pbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的5bit SI 2次置换
            int[][] S5table = ReadShiftInvariantTT(S5binfile, 5);
            //读取所有有直接TI的3bit SI 2次置换
            int[][] S3table = ReadShiftInvariantTT(S3binfile, 3);
            //读取所有有直接TI的8bit SI 线性置换
            int[][] Ptable = ReadShiftInvariantTT(Pbinfile, 8);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            long length = S5table.Length * S3table.Length * Ptable.Length;
            for (long num = S5table.Length * S3table.Length; num < length; num++)
            {
                if ((num & 0xff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0 * num / length);
                }
                //从num中获取P的序号
                int Pno = (int)(num / (S5table.Length * S3table.Length));
                int Sno = (int)(num % (S5table.Length * S3table.Length));
                int Sno5 = Sno / S3table.Length;
                int Sno3 = Sno % S3table.Length;
                //从P中获取置换表
                //int[] Ptable = GetBitPerm8(P);

                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_ShiftInvariant_53SPN(table, Ptable[Pno], S5table[Sno5], S3table[Sno3]);
                }
                //if (!CheckPermutaion(table))
                //{
                //    System.Console.WriteLine("None permuataion!");
                //    return errorTable;
                //};
                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    mindiff = diff;

                if ((num & 0xffff) == 0)
                    System.Console.WriteLine("diff={0},num={1:x}", diff, num);

                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    //System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear > maxnL)
                        maxnL = nonLinear;
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数={0:x}", num);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();

                    }

                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();


            return errorTable;
        }
        //检查P是否会将两个Sbox部分交互
        public bool Check53P(int[] Ptable)
        {
            for (int i = 0; i < 3; i++)
            {
                if (Ptable[i] >= 3)
                    return true;
            }
            return false;
        }

        //5bit/3bit SI Sbox+8bit BitP 线性变换
        public int[] SearchOptimal_ShiftInvariant_53SBitP(string S5binfile, string S3binfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的5bit SI 2次置换
            int[][] S5table = ReadShiftInvariantTT(S5binfile, 5);
            //读取所有有直接TI的3bit SI 2次置换
            int[][] S3table = ReadShiftInvariantTT(S3binfile, 3);



            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            long length = S5table.Length * S3table.Length * 40320;
            for (long num = S5table.Length * S3table.Length; num < length; num++)
            {
                if ((num & 0xffff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0 * num / length);
                }
                //从num中获取P的序号
                int Pno = (int)(num / (S5table.Length * S3table.Length));
                int Sno = (int)(num % (S5table.Length * S3table.Length));
                int Sno5 = Sno / S3table.Length;
                int Sno3 = Sno % S3table.Length;
                //从P中获取置换表
                int[] Ptable = GetBitPerm8(Pno);
                //检查是否P会交互两个S
                if (!Check53P(Ptable))
                    continue;

                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_ShiftInvariant_53SBitP(table, Ptable, S5table[Sno5], S3table[Sno3]);
                }
                //if (!CheckPermutaion(table))
                //{
                //    System.Console.WriteLine("None permuataion!");
                //    return errorTable;
                //};
                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    mindiff = diff;

                if ((num & 0xffff) == 0)
                    System.Console.WriteLine("mindiff={0},num={1:x}", mindiff, num);

                if (diff <= optimalDiff)
                {
                    int nonLinear = Nonlinear(table, size);
                    //System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                    if (nonLinear > maxnL)
                        maxnL = nonLinear;
                    if (nonLinear >= optimalNonlinear)
                    {
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < 256; j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数={0:x}", num);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();

                    }

                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();


            return errorTable;
        }

        //由bitP构造P的table
        public int[] GetBitPTable(int[] P)
        {
            int[] Ptable = new int[0x1 << P.Length];
            for (int i = 0; i < Ptable.Length; i++)
            {
                Ptable[i] = 0;
                for (int j = 0; j < 8; j++)
                {
                    if (GetBit(i, P[j]) == 1)
                      Ptable[i] = SetBit(Ptable[i], j);
                }
            }
            return Ptable;
        }

        //5bit/3bit SI Sbox+8bit BitP 线性变换+固定常数
        public int[] SearchOptimal_ShiftInvariant_53SBitP_FixC(string S5binfile, string S3binfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的5bit SI 2次置换
            int[][] S5table = ReadShiftInvariantTT(S5binfile, 5);
            //读取所有有直接TI的3bit SI 2次置换
            int[][] S3table = ReadShiftInvariantTT(S3binfile, 3);
            //int mindiff5 = 32;
            //for (int i = 0; i < S5table.Length; i++)
            //{
            //    int diff = DiffUniform(S5table[i], 5);
            //    System.Console.WriteLine("{1}:\t{0}\t{2}", diff,i);
            //    if (diff < mindiff5)
            //        mindiff5 = diff;
            //}
            //System.Console.WriteLine("Sbox5:\tmindiff={0}", mindiff5);
            //int mindiff3 = 8;
            //for (int i = 0; i < S3table.Length; i++)
            //{
            //    int diff = DiffUniform(S3table[i], 3);
            //    System.Console.WriteLine("{1}:\t{0}", diff, i);
            //    if (diff < mindiff3)
            //        mindiff3 = diff;
            //}
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);

            long length = 40320 * (S5table.Length * S3table.Length);
            for (long num = 0; num < length; num++)
            {
                if ((num & 0xffff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%,mindiff={2}", num, 100.0 * num / length,mindiff);
                }
                //从num中获取P的序号
                //int Pno = (int)(num / (S5table.Length * S3table.Length));
                //int Sno = (int)(num % (S5table.Length * S3table.Length));
                //int Sno5 = Sno / (S3table.Length);
                //int Sno3 = Sno % (S3table.Length);

                int Sno = (int)(num / 40320);
                int Pno = (int)(num % 40320);
                int Sno5 = Sno / (S3table.Length);
                int Sno3 = Sno % (S3table.Length);
              


                //从P中获取置换表
                int[] Ptable = GetBitPerm8(Pno);
                //检查是否P会交互两个S
                if (!Check53P(Ptable))
                    continue;
                //将bit置换做表
                Ptable = GetBitPTable(Ptable);
                for (int C = 1; C < 2; C++)
                {
                    //table初始化为恒等变换
                    for (int i = 0; i < len; i++)
                    {
                        table[i] = i;
                    }

                    for (int i = 0; i < round; i++)
                    {
                        OneRoundTrans_ShiftInvariant_53SBitP_FixC(table, Ptable, S5table[Sno5], S3table[Sno3], C);
                    }
                    //if (!CheckPermutaion(table))
                    //{
                    //    System.Console.WriteLine("None permuataion!");
                    //    return errorTable;
                    //};
                    //求解最终S盒的差分均匀度、非线性度
                    int diff = DiffUniform(table, size);
                    //sw.WriteLine("{0}\t{1}",Pno,diff);
                   
                    if (diff < mindiff)
                        mindiff = diff;

                    if(C==1 &&((num&0xff)==0))
                        System.Console.WriteLine("num={0:x},C={3},percent={1}%,mindiff={2}", num, 100.0 * num / length, mindiff,C);

                    if (diff <= optimalDiff)
                    {
                        int nonLinear = Nonlinear(table, size);
                        //System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        if (nonLinear > maxnL)
                            maxnL = nonLinear;
                        if (nonLinear >= optimalNonlinear)
                        {
                            System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                            count++;
                            sw.WriteLine("******************");
                            sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                            sw.WriteLine("输出表：");
                            for (int j = 0; j < 256; j++)
                                sw.Write("{0:x2}\t", table[j]);
                            sw.WriteLine("\n");
                            sw.WriteLine("内部函数={0:x}", num);
                            sw.WriteLine("******************");
                            sw.WriteLine("\n");
                            sw.WriteLine("\n");
                            sw.WriteLine("\n");
                            sw.Flush();

                        }

                    }
                    //差分大于16的就不考虑常数了
                    if (C == 0 && diff > 16)
                        break;
                }
            }
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();


            return errorTable;
        }

        public int RotatedShift(int x, int s, int bitlen)
        {
            int sl = s % bitlen;
            int mask=(0x1<<bitlen)-1;
            return ((x << sl) & mask) | (x >> (bitlen - sl));
        }
            //5bit/3bit SI Sbox+8bit BitP 线性变换+固定常数
        public void SearchOptimal_ShiftInvariant_53SBitP_ShiftC(string S5binfile, string S3binfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            
            int count = 0;
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的5bit SI 2次置换
            int[][] S5table = ReadShiftInvariantTT(S5binfile, 5);
            //读取所有有直接TI的3bit SI 2次置换
            int[][] S3table = ReadShiftInvariantTT(S3binfile, 3);

            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            long length = 40320 * (S5table.Length * S3table.Length);
            Parallel.For(0, length, new ParallelOptions { MaxDegreeOfParallelism = 4 }, num =>
                {
                    if ((num & 0xffff) == 0)
                    {
                        System.Console.WriteLine("num={0:x},percent={1}%,mindiff={2}", num, 100.0 * num / length, mindiff);
                    }
                    int Sno = (int)(num / 40320);
                    int Pno = (int)(num % 40320);
                    int Sno5 = Sno / (S3table.Length);
                    int Sno3 = Sno % (S3table.Length);
                    //从P中获取置换表
                    int[] Ptable = GetBitPerm8(Pno);
                    //检查是否P会交互两个S
                    if (!Check53P(Ptable))
                        return;
                    //将bit置换做表
                    Ptable = GetBitPTable(Ptable);
                    int[] table = new int[len];
                    //table初始化为恒等变换
                    for (int i = 0; i < len; i++)
                    {
                        table[i] = i;
                    }
                    for (int i = 0; i < round; i++)
                    {
                        OneRoundTrans_ShiftInvariant_53SBitP_FixC(table, Ptable, S5table[Sno5], S3table[Sno3], RotatedShift(0x78,i, size));
                    }
                    //求解最终S盒的差分均匀度、非线性度
                    int diff = DiffUniform(table, size);
                    //sw.WriteLine("{0}\t{1}",Pno,diff);

                    if (diff < mindiff)
                        Interlocked.Exchange(ref mindiff, diff);


                    if (diff <= optimalDiff)
                    {
                        int nonLinear = Nonlinear(table, size);
                        //System.Console.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        if (nonLinear > maxnL)
                            Interlocked.Exchange(ref maxnL, nonLinear); 
                        if (nonLinear >= optimalNonlinear)
                        {
                            lock(sw)
                            {
                            System.Console.WriteLine("*********diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                            count++;
                            sw.WriteLine("******************");
                            sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                            sw.WriteLine("输出表：");
                            for (int j = 0; j < 256; j++)
                                sw.Write("{0:x2}\t", table[j]);
                            sw.WriteLine("\n");
                            sw.WriteLine("内部函数={0:x}", num);
                            sw.WriteLine("******************");
                            sw.WriteLine("\n");
                            sw.WriteLine("\n");
                            sw.WriteLine("\n");
                            sw.Flush();
                            }
                        }

                    }
                });

            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();


            return;
        }
    }
}
