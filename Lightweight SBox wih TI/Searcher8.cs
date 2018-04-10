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
        //Basic Tools
        #region
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
                //置bit操作
        void SetBit(byte[] input, int pos)
        {
            input[pos/8]=(byte)(input[pos/8]|(0x1<<(pos%8)));
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

        //计算HW
        public int HW(int x, int len)
        {
            int c = 0;
            for(int i=0;i<len;i++)
                if(GetBit(x,i)!=0)
                   c++;
            return c;
        }
        public int RotatedShift(int x, int s, int bitlen)
        {
            int sl = s % bitlen;
            int mask = (0x1 << bitlen) - 1;
            return ((x << sl) & mask) | (x >> (bitlen - sl));
        }
        #endregion

        //Crypto Tools
        #region
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
        public int Nonlinear(int[] sbox, int size)
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
        public int[] MoebiusTrans(int varNum, byte[] func)
        {
            int[] ANF = new int[0x1 << varNum];
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
            return ANF;
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
        public int[] InvMoebiusTrans1(int varNum, int[] ANF)
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

            return ANF1;
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
        //根据ANF算最大次数,v表示变量个数
        public int MaxDegree(int[] ANF, int v)
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
        public int[] CreatePermTable(int l, int[] Perm)
        {
            int len = (int)Math.Pow(2, l);
            int[] table = new int[len];
            for (int i = 0; i < len; i++)
            {
                table[i] = 0;
                for (int j = 0; j < l; j++)
                    if (GetBit(i, j) != 0)
                        table[i] = SetBit(table[i], Perm[j]);
            }

            //检查
            //Array.Sort(table);
            //for (int i = 0; i < len - 1; i++)
            //    if (table[i] == table[i + 1])
            //        return null;
            //检查结束
            return table;
        }

        //Linear Transformation
        public int LinearP(int[][] PMatrix, int x,int size)
        {
            int y=0;
            for (int i = 0; i < size; i++)
            {
                int t = 0;
                for (int j = 0; j < size; j++)
                {
                    t = PMatrix[i][j] * GetBit(x, j) ^ t;
                }
                if (t == 1)
                    y=SetBit(y, i);
            }
            return y;
        }
        //Linear Transformation:Field Mulitiplication
        public int FieldMultiplication(int x,int FieldNo,int MultiNo)
        {
            int y = 0;
            for (int i = size - 1; i > -1; i--)
            {
                if (GetBit(MultiNo, i) > 0)
                    y = y ^ x;
                if(i!=0)
                    y = FieldM2(y, FieldNo);
            }
            return y;
        }
        //Field M2
        public int FieldM2(int x, int FieldNo)
        {
            int y = 0;
            if ((x >> (size - 1)) > 0)
                y = ((x << 1) & ((0x1 << size) - 1)) ^ ((FieldNo << 1) ^ 0x01);
            else
                y = ((x << 1) & ((0x1 << size) - 1));
            return y;
        }
        //Transform Truth Table from int[] to bitwise byte[]
        public byte[] TTTransformation(int[] TT)
        {
            byte[] TTnew=new byte[TT.Length/8];
            for(int i=0;i<TT.Length;i++)
            {
                if(TT[i]==1)
                    SetBit(TTnew,i);
            }
            return TTnew;
        }
        //Transform Truth Table from int[] to bitwise byte[]
        public int[] TTTransformation1(byte[] TT)
        {
            int[] TTnew = new int[TT.Length* 8];
            for (int i = 0; i < TTnew.Length; i++)
            {
                if (GetBit(TT, i)== 1)
                    TTnew[i]=1;
            }
            return TTnew;
        }
        #endregion

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
        //二元阵求逆
        #region
        //
        //函数GJ_elimination，对0-1矩阵进行高斯消元
        public void GJ_elimination(int[][] M, int rownum, int colnum)
        {
            int i, j, k, flag, temp;

            int rank = 0;

            int startcol, currrow, postrow = 0;

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
                tempM[i] = new int[dim * 2];
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

        //Specific tools
        #region
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

        //根据SI规则构造完整置换表
        public int[] ShiftInvariantTT(byte[] TT, int bitlen)
        {
            int len = 0x1 << bitlen;
            int mask = len - 1;
            int[] table=new int[len];
            for (int x = 0; x < len; x++)
            {
                table[x]=0;
                int xt = x;
                for (int bit = 0; bit < bitlen; bit++)
                {
                    if (GetBit(TT, xt) == 1)
                          table[x] = SetBit(table[x], bit); 
                    //table[x] = (table[x] ^ 0x1);
                    //if (bit != bitlen - 1)
                    //    table[x] = table[x] << 1;

                    // xt = ((xt << 1)&mask) | (xt>>(bitlen-1));
                    xt = ((xt << (bitlen - 1)) & mask) | (xt >> 1);
                }
            }
            return table;
        }

        //从SI TI的bit真值表中构造置换表
        public int[][] ReadShiftInvariantTT_Sbox(string TTfile, int bitlen)
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
        //读取SI TI的bit真值表
        public byte[][] ReadShiftInvariantTT_OneBit(string TTfile, int bitlen)
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
            byte[][] table = new byte[num][];
            //读写指针还原
            br.BaseStream.Seek(0, SeekOrigin.Begin);
            for (int i = 0; i < num; i++)
            {
                table[i] = br.ReadBytes(bytelen);
            }
            br.Close();
            fs.Close();
            return table;
        }
        //根据编号构造*2的线性变换
        public int[][] BuildM2PMatrix(int Pno, int size)
        {
            int[][] Pmatrix = new int[size][];
            for (int i = 0; i < size; i++)
            {
                Pmatrix[i] = new int[size];
                for (int j = 0; j < size; j++)
                    if (j == (i + 1) % size)
                        Pmatrix[i][j] = 1;
            }
            //根据Pno确定第一列前n-1个值
            for (int i = 0; i < size - 1; i++)
            {
                Pmatrix[i][0] = GetBit(Pno, i);
            }
            return Pmatrix;
        }
        //根据编号构造任意线性变换
        public int[][] BuildAnyPMatrix(long Pno, int size)
        {
            int[][] Pmatrix = new int[size][];
            for (int i = 0; i < size; i++)
            {
                Pmatrix[i] = new int[size];
                for (int j = 0; j < size; j++)
                        Pmatrix[i][j] = GetBit(Pno,i*size+j);
            }
            return Pmatrix;
        }

        //根据ANF评估代价
        //ANF中AND的数量和XOR的数量
        //cost[0]:AND的数量;cost[1]:XOR的数量
        public int[] ANFCost_Hardware(byte[] TT,int size)
        {
            int[] ANF = MoebiusTrans(size, TT);
            int[] Cost = new int[2];
            Cost[0] = 0;
            Cost[1] = -1;
            for (int i = 0; i < ANF.Length; i++)
            {
                if (ANF[i] == 1)
                {
                    Cost[1]++;
                    if (HW(i,size) > 1)
                        Cost[0]++;
                }
            }
            Cost[0] = Cost[0] * size;
            Cost[1] = Cost[1] * size;
            return Cost;
        }
        public int[] ANFCost_Hardware_1b(byte[] TT, int size)
        {
            int[] ANF = MoebiusTrans(size, TT);
            int[] Cost = new int[2];
            Cost[0] = 0;
            Cost[1] = -1;
            for (int i = 0; i < ANF.Length; i++)
            {
                if (ANF[i] == 1)
                {
                    Cost[1]++;
                    if (HW(i, size) > 1)
                        Cost[0]++;
                }
            }
           
            return Cost;
        }
        //根据ANF评估代价
        //ANF中AND的数量和XOR的数量
        //cost[0]:AND的数量;cost[1]:XOR的数量;cost[2]:Shift的数量；cost[3]:其他指令
        public int[] ANFCost_Software(byte[] TT,int size)
        {
            int[] ANF = MoebiusTrans(size, TT);
            int[] Cost = new int[4];
            Cost[0] = 0;
            Cost[1] = -1;
            Cost[2] = 0;
            Cost[3] = 0;
            for (int i = 0; i < ANF.Length; i++)
            {
                if (ANF[i] == 1)
                {
                    Cost[1]++;
                    if (HW(i, size) > 1)//2次项
                    {
                        Cost[0]++;
                        if ((i & 0x01) == 0)//没有x0
                            Cost[2] = Cost[2] + 2;
                        else//有一个是x0
                            Cost[2] = Cost[2] + 1;
                    }
                    if (HW(i, size) == 1)//1次项
                    {
                        if ((i & 0x01) == 0)//不是x0
                            Cost[2] = Cost[2] + 1;
                    }
                }
            }
            return Cost;
        }
        //根据多项式评估硬件*4的代价
        //返回所用的XOR数量
        public int M4Cost_Hardware(int Pno)
        {
            return 2 * HW(Pno,size);
        }
        //根据多项式评估软件*4的代价
        //返回所用的XOR数量
        public int[] M4Cost_Software(int Pno)
        {
            int[] Cost = new int[4];
            Cost[0] = 0;
            Cost[1] = 2;//1个Xor
            Cost[2] = 2;//1个循环移位
            Cost[3] = 4;//Branch的1个TEST和一个JMP
            return Cost;
        }
        //评估整体unprotecked Sbox的实现代价
        public void WriteSIM4PCost_Unprotected(StreamWriter sw,byte[] TT,int Pno)
        {
            sw.WriteLine("Worst Case Unprotected Raw Implementation:");
            int[] HardwareCost=ANFCost_Hardware(TT,size);
            int[] SoftwareCost=ANFCost_Software(TT,size);
            int[] Temp=M4Cost_Software(Pno);
            for(int i=0;i<4;i++)
                SoftwareCost[i]+=Temp[i];
            double sumHard=HardwareCost[0]*1.67+HardwareCost[1]*2.67;
            double sumSoft_Rs=SoftwareCost[0]+SoftwareCost[1]+SoftwareCost[2]+SoftwareCost[3];
            double sumSoft_NoRs=SoftwareCost[0]+SoftwareCost[1]+SoftwareCost[2]*3+SoftwareCost[3];
            HardwareCost[1]=HardwareCost[1]+M4Cost_Hardware(Pno);
            sw.WriteLine("Hardware(In theory): AND={0},XOR={1}",HardwareCost[0]*round,HardwareCost[1]*round);
            sw.WriteLine("Overall GE={0}", sumHard * round);
            sw.WriteLine("Software Cycles: AND={0},XOR={1},RotatedShift={2},Other={3}", SoftwareCost[0] * round, SoftwareCost[1] * round, SoftwareCost[2] * round, SoftwareCost[3] * round);
            sw.WriteLine("Overall cycles with Rotated Shift={0}, without Rotated Shift={1}", sumSoft_Rs * round, sumSoft_NoRs * round);
            sw.WriteLine("Worst Case Unprotected Round-based Implementation:");
            sw.WriteLine("Hardware(In theory): AND={0},XOR={1}", HardwareCost[0], HardwareCost[1]);
            sw.WriteLine("Overall GE={0}", sumHard);
            sw.WriteLine("Software Cycles: AND={0},XOR={1},RotatedShift={2},Other={3}", SoftwareCost[0], SoftwareCost[1], SoftwareCost[2], SoftwareCost[3]);
            sw.WriteLine("Overall cycles with Rotated Shift={0}, without Rotated Shift={1}", sumSoft_Rs, sumSoft_NoRs );
            sw.Flush();
        }
        //评估整体unprotecked Sbox的实现代价
        public void WriteSIM4PCost_TI1b(StreamWriter sw, byte[] TT1b, int Pno)
        {
            int shares = 3;
            sw.WriteLine("Worst Case TI Raw Implementation:");
            int[] HardwareCost = ANFCost_Hardware(TT1b,size*(shares-1));
            HardwareCost[1] = HardwareCost[1] + shares * M4Cost_Hardware(Pno);
            int[] SoftwareCost = ANFCost_Software(TT1b, size * (shares - 1));
            int[] Temp = M4Cost_Software(Pno);
            for (int i = 0; i < 4; i++)
                SoftwareCost[i] += Temp[i];
            double sumHard = HardwareCost[0] * 1.67 + HardwareCost[1] * 2.67;
            double sumSoft_Rs = SoftwareCost[0] + SoftwareCost[1] + SoftwareCost[2] + SoftwareCost[3];
            double sumSoft_NoRs = SoftwareCost[0] + SoftwareCost[1] + SoftwareCost[2] * 3 + SoftwareCost[3];
            
            sw.WriteLine("Hardware(In theory): AND={0},XOR={1}", HardwareCost[0] * round, HardwareCost[1] * round);
            sw.WriteLine("Overall GE={0}", sumHard * round);
            sw.WriteLine("Software Cycles: AND={0},XOR={1},RotatedShift={2},Other={3}", SoftwareCost[0] * round, SoftwareCost[1] * round, SoftwareCost[2] * round, SoftwareCost[3] * round);
            sw.WriteLine("Overall cycles with Rotated Shift={0}, without Rotated Shift={1}", sumSoft_Rs * round, sumSoft_NoRs * round);
            
            sw.WriteLine("Worst Case TI 1 bit Implementation:");
            HardwareCost = ANFCost_Hardware_1b(TT1b, size * (shares - 1));
            HardwareCost[1] = HardwareCost[1] + shares * M4Cost_Hardware(Pno);
            sumHard = HardwareCost[0] * 1.67 + HardwareCost[1] * 2.67;
            sw.WriteLine("Hardware(In theory): AND={0},XOR={1}", HardwareCost[0], HardwareCost[1]);
            sw.WriteLine("Overall GE={0}", sumHard);
            sw.WriteLine("Software Cycles: AND={0},XOR={1},RotatedShift={2},Other={3}", SoftwareCost[0], SoftwareCost[1], SoftwareCost[2], SoftwareCost[3]);
            sw.WriteLine("Overall cycles with Rotated Shift={0}, without Rotated Shift={1}", sumSoft_Rs, sumSoft_NoRs);
            sw.Flush();
        }
        //用Verilog写出未保护的Sbox
        public void Print_FullSbox_Unprotected(string path, int[] Sbox,long num,string sname)
        {
            String filename=String.Format( path+sname+"{0}_R{1}_{2}.v",size,round,num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            
            sw.WriteLine("module S{0}_R{1}_{2}(in,out);",size,round,num );
            sw.WriteLine("  input[{0}:0] in;", size - 1);
            sw.WriteLine("  output[{0}:0] out;", size - 1);
            sw.WriteLine("  reg[{0}:0] out;", size - 1);
            sw.WriteLine("  always@(in)");
            sw.WriteLine("  begin");
            sw.WriteLine("    case(in)");
            for (int i = 0; i < Sbox.Length; i++)
                sw.WriteLine("   {0}: out<={1};", i, Sbox[i]);
            sw.WriteLine(" endcase");
            sw.WriteLine("end");
            sw.WriteLine("endmodule");
            sw.Close();
            fs.Close();
        }

        //用Verilog写出未保护的Sbox的一轮
        public void Print_OneRoundSbox_Unprotected(string path, int[] Sbox, long num, string sname)
        {
            String filename = String.Format(path + sname + "{0}_R{1}_{2}.v", size, round, num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);

            sw.WriteLine("module S{0}_R{1}_{2}(in,out);", size, round, num);
            sw.WriteLine("  input[{0}:0] in;", size - 1);
            sw.WriteLine("  output[{0}:0] out;", size - 1);
            sw.WriteLine("  reg[{0}:0] out;", size - 1);
            sw.WriteLine("  always@(in)");
            sw.WriteLine("  begin");
            sw.WriteLine("    case(in)");
            for (int i = 0; i < Sbox.Length; i++)
                sw.WriteLine("   {0}: out<={1};", i, Sbox[i]);
            sw.WriteLine(" endcase");
            sw.WriteLine("end");
            sw.WriteLine("endmodule");
            sw.Close();
            fs.Close();
        }

        //用Verilog写出1bit SI
        public void Print_SI1bit_Unprotected(string path, int[] Sbox, long num, string sname)
        {
            String filename = String.Format(path + sname + "{0}_R{1}_{2}.v", size, round, num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);

            sw.WriteLine("module S{0}_R{1}_{2}(in,out);", size, round, num);
            sw.WriteLine("  input[{0}:0] in;", size - 1);
            sw.WriteLine("  output out;");
            sw.WriteLine("  reg out;");
            sw.WriteLine("  always@(in)");
            sw.WriteLine("  begin");
            sw.WriteLine("    case(in)");
            for (int i = 0; i < Sbox.Length; i++)
                sw.WriteLine("   {0}: out<={1};", i, Sbox[i]&0x01);
            sw.WriteLine(" endcase");
            sw.WriteLine("end");
            sw.WriteLine("endmodule");
            sw.Close();
            fs.Close();
        }
        //写出特定Sbox对应的脚本

        //用Verilog写出1bit SI
        public void Print_TI1b(string path, int[] Sbox, long num, string sname,int shares)
        {
            String filename = String.Format(path + sname + "{0}_R{1}_{2}.v", size, round, num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);

            sw.WriteLine("module STI{0}_R{1}_{2}(in,out);", size, round, num);
            sw.WriteLine("  input[{0}:0] in;", (shares-1)*size - 1);
            sw.WriteLine("  output out;");
            sw.WriteLine("  reg out;");
            sw.WriteLine("  always@(in)");
            sw.WriteLine("  begin");
            sw.WriteLine("    case(in)");
            for (int i = 0; i < Sbox.Length; i++)
                sw.WriteLine("   {0}: out<={1};", i, Sbox[i] & 0x01);
            sw.WriteLine(" endcase");
            sw.WriteLine("end");
            sw.WriteLine("endmodule");
            sw.Close();
            fs.Close();
        }
 
        //Convert ANF to ANF terms
        public ANFterm[] ConvertANF(int[] ANF, int shares)
        {
            
            int count = 0;
            for (int i = 0; i < ANF.Length; i++)
            {
                if (ANF[i] == 1)
                {
                    count++;
                }
            }
            ANFterm[] res = new ANFterm[count];
            int n = 0;
            for (int i = 0; i <ANF.Length; i++)
            {
                if (ANF[i] == 1)
                {
                    int d=0;
                    //Get the bit index
                    for (int j = 0; j < size*(shares-1); j++)
                    {
                        if (GetBit(i, j) == 1)
                            d++;
                    }
                    //Allocate space
                    int[] ind=new int[d];
                    int m=0;
                    for (int j = 0; j < size * (shares - 1); j++)
                    {
                        if (GetBit(i, j) == 1)
                        {
                           // int bit = j / (shares - 1);
                          //  ind[m]=j+bit*2;
                            ind[m] = j;
                            m++;
                        }
                    }
                    //Add term
                    res[n] = new ANFterm(ind);
                    n++;
                }
            }
            return res;
        }
        //Write the S function in Sbox
        public void WriteASM_S(StreamWriter sw, int shares, byte[] TT)
        {
            int[] ANF = MoebiusTrans(size * (shares - 1), TT);
            //Compute ANF terms
            ANFterm[] anf = ConvertANF(ANF, shares);

            int len_r = size * (shares - 1);
            int len_e = size * (shares + 1);
            sw.WriteLine(".func TI_S");
            sw.WriteLine("TI_S:");
            sw.WriteLine("push {lr}");
            sw.WriteLine("push {r4-r7}");
            for (int i = 0; i < 8; i++)
                sw.WriteLine("nop");
            //r0 is the input register, contains shared data like x00x01x02x00x10x11x12x10....
            //Compute the shift version
            //Pass r0 as input
            //Pass r1 as ShareMask2
            //Pass r2 as ShareMask1
            //Pass r3 as AndMask
            //Note that and/eor/lsl/lsr do not work with high regs

            //Use r{0} as temp
            //Use r{1} as temp1
            //Use r{2} as temp2
            //Use r{3-(len_r+2)} as shift results 0-7
            //Use r{len_r+3} as AndMask:0x0000ffff
            //Can not use r13 in Cortex M0, so leave r{len_r+5}
            //Use r{len_r+6} as ShareMask2:0b'1100110011001100
            //Use r{len_r+4} as ShareMask1:0b'0010001000100010
            //Save Masks
            sw.WriteLine("   mov r{0},r3", len_r + 3);
            sw.WriteLine("   mov r{0},r1", len_r + 6);
            sw.WriteLine("   mov r{0},r2", len_r + 4);
            //Get Shift results
            sw.WriteLine("   mov r3,r0");
            for (int i = 1; i < len_r; i++)
            {
                if ((i % (shares - 1) == 0) && (i > 1))
                {
                    //Clean R0 to avoid leakage
                    sw.WriteLine("   mov r0,#0");//
                    //Right shift varables=Rotated shift by 4
                    sw.WriteLine("   mov r0,r{0}", i + 1);//temp=reg[i+1];
                    sw.WriteLine("   lsl r1,r0,#{0}", len_e - 4);//temp1=temp<<12
                    sw.WriteLine("   lsr r0,#4");//temp=temp>>4;
                    sw.WriteLine("   eor r0,r1");//temp=temp^temp1;
                    if (len_e != 32)
                    {
                        sw.WriteLine("   mov r1,r{0}", len_r + 3);//temp1=Andmask;
                        sw.WriteLine("   and r0,r1");//temp=temp&Andmask;
                    }
                    sw.WriteLine("   mov r{0},r0", i + 3);//r{i+3}=temp;
                }
                else
                {
                    //Right shift shares by 1=Rotated shift shares by 1
                    sw.WriteLine("   mov r0,r{0}", i + 2);//temp=reg[i+2];
                    //Get high two bits
                    sw.WriteLine("   mov r1,r{0}", len_r + 6);//temp1=ShareMask2;
                    sw.WriteLine("   and r1,r0");//temp1=reg[i+2]&ShareMask2
                    sw.WriteLine("   lsr r1,#1");//temp1=temp1>>1
                    //Change to low 1 bit
                    sw.WriteLine("   mov r2,r{0}", len_r + 4);//temp2=ShareMask1;
                    sw.WriteLine("   mov r0,r{0}", i + 2);//temp=reg[i+2];
                    sw.WriteLine("   and r0,r2");//temp=ShareMask1&reg[i+2];
                    sw.WriteLine("   lsr r0,#1");//temp=temp>>1
                    //Left shift 3
                    sw.WriteLine("   lsl r2,r0,#3");//temp2=temp<<3;
                    sw.WriteLine("   eor r0,r2");//temp=temp^(temp<<3);
                    //Paste together
                    sw.WriteLine("   eor r0,r1");//temp=temp2^temp1
                    sw.WriteLine("   mov r{0},r0", i + 3);//r{i+3}=temp;
                }
            }
            //Sort anf as degree and delta_index
            Array.Sort(anf);
            //Use r{0} as result
            //Use r{1} as temp1
            //Use r{2} as temp2
            sw.WriteLine("   mov r0,#0");
            for (int i = 0; i < anf.Length; i++)
            {
                if (anf[i].d == 1)//linear term
                {
                    sw.WriteLine("   mov r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    sw.WriteLine("   eor r0,r1");
                }
                else//only deal with degree 2 funcitons
                {
                    sw.WriteLine("   mov r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    sw.WriteLine("   mov r2,r{0}", anf[i].base_index + anf[i].other_index[0] + 3);//temp2=reg[];
                    sw.WriteLine("   and r1,r2");//temp1=temp1&temp2
                    //Xor temp to result
                    sw.WriteLine("   eor r0,r1");

                }
            }
            for (int i = 0; i < 8; i++)
                sw.WriteLine("nop");
            sw.WriteLine("pop {r4-r7}");
            sw.WriteLine("pop {pc}");
            sw.WriteLine(".endfunc");
        }
        //Write the S function in Sbox
        public void WriteASM_S_SharesGroup(StreamWriter sw, int shares, byte[] TT)
        {
            int[] ANF = MoebiusTrans(size * (shares - 1), TT);
            //Compute ANF terms
            ANFterm[] anf = ConvertANF(ANF, shares);

            int len_r = size * (shares - 1);
            int len_e = 2*size * (shares);
            sw.WriteLine(".func TI_S");
            sw.WriteLine("TI_S:");
            sw.WriteLine("push {lr}");
            sw.WriteLine("push {r4-r7}");
            for (int i = 0; i < 8; i++)
                sw.WriteLine("nop");
            //r0 is the input register, contains shared data like x00x10x20x30....
            //Compute the shift version
            //Pass r0 as input
            //Pass r1 as GetMask
            //Pass r2 as AndMask
            //Note that and/eor/lsl/lsr do not work with high regs

            //Use r{0} as temp
            //Use r{1} as temp1
            //Use r{2} as temp2
            //Use r{3-(len_r+2)} as shift results 0-7
            //Use r{len_r+3} as AndMask:0x00ffffff
            //Can not use r13 in Cortex M0, so leave r{len_r+5}
            //Use r{len_r+4} as GetMask:0x000f0f0f
            //Save Masks
            sw.WriteLine("   mov r{0},r1", len_r + 4);
            sw.WriteLine("   mov r{0},r2", len_r + 3);
            //Get Shift results
            sw.WriteLine("   mov r3,r0");
            for (int i = 1; i < len_r; i++)
            {
                if (i % (size) == 0)//shift shares
                {
                    //Right shift shares==shift 2*size
                    sw.WriteLine("   mov r0,r{0}", i - size + 3);//temp=reg[i -size+3];
                    sw.WriteLine("   lsl r1,r0,#{0}", len_e - 2 * size);//temp1=temp<<16
                    sw.WriteLine("   lsr r0,#{0}",2*size);//temp=temp>>8;
                    sw.WriteLine("   eor r0,r1");//temp=temp^temp1;
                    if (len_e != 32)
                    {
                        sw.WriteLine("   mov r1,r{0}", len_r + 3);//temp1=Andmask;
                        sw.WriteLine("   and r0,r1");//temp=temp&Andmask;
                    }
                    sw.WriteLine("   mov r{0},r0", i + 3);//r{i+3}=temp;

                }
                else//shift bits
                {
                    ///Right shift bits by 1
                    sw.WriteLine("   mov r0,r{0}", i + 2);//temp=reg[i+2];
                    sw.WriteLine("   lsr r1,r0,#1");//temp1=temp>>1
                    sw.WriteLine("   mov r2,r{0}", len_r + 4);//temp2=GetMask;
                    sw.WriteLine("   and r2,r1");//temp2=temp1&GetMask
                    sw.WriteLine("   lsl r0,r2,#{0}",size);//temp=temp2<<4
                    sw.WriteLine("   eor r0,r2");//temp=temp^temp2
                    sw.WriteLine("   mov r{0},r0", i + 3);//r{i+3}=temp;
                }
            }
            //Sort anf as degree and delta_index
            Array.Sort(anf);
            //Use r{0} as result
            //Use r{1} as temp1
            //Use r{2} as temp2
            sw.WriteLine("   mov r0,#0");
            for (int i = 0; i < anf.Length; i++)
            {
                if (anf[i].d == 1)//linear term
                {
                    sw.WriteLine("   mov r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    sw.WriteLine("   eor r0,r1");
                }
                else//only deal with degree 2 funcitons
                {
                    sw.WriteLine("   mov r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    sw.WriteLine("   mov r2,r{0}", anf[i].base_index + anf[i].other_index[0] + 3);//temp2=reg[];
                    sw.WriteLine("   and r1,r2");//temp1=temp1&temp2
                    //Xor temp to result
                    sw.WriteLine("   eor r0,r1");

                }
            }
            for (int i = 0; i < 8; i++)
                sw.WriteLine("nop");
            sw.WriteLine("pop {r4-r7}");
            sw.WriteLine("pop {pc}");
            sw.WriteLine(".endfunc");
        }
        //Write the P function in Sbox
        public void WriteASM_P(StreamWriter sw,int shares)
        {
            sw.WriteLine(".func TI_P");
            sw.WriteLine("TI_P:");
            sw.WriteLine("push {lr}");
            sw.WriteLine("push {r4-r5}");
            int len_e = size * (shares + 1);
            //r0 as the sharedinput
            //r1 as the address of PMaskTable
            //r2 as AndMask
            //r3 as temp1
            //r4 as temp2
            //Do Pmatrix
            sw.WriteLine("   lsl r3,r0,#4");//temp1=result<<4;
            sw.WriteLine("   and r3,r2 ");//temp1=temp1&AndMask;
            sw.WriteLine("   lsr r0,#{0}",len_e - 4);//result=result>>12;
            sw.WriteLine("   mov r4,#7");//cannot and a constant, use a reg instead
            sw.WriteLine("   and r0,r4");//result=result&0x07;
            sw.WriteLine("   lsl r0,#2");//result=result*4;
            sw.WriteLine("   ldr r0,[r1,r0]");//result=PMTable[resullt];
            sw.WriteLine("   eor r0,r3");//result=temp1^result;

            sw.WriteLine("pop {r4-r5}");
            sw.WriteLine("pop {pc}");
            sw.WriteLine(".endfunc");
        }

        //用ASM写出软件实现TI的代码
        public void Print_TI1b_ASM(string path, long num, string sname, int shares,byte[] TT)
        {
            
            //Write the ASM
            String filename = String.Format(path + sname + "{0}_R{1}_{2}.S", size, round, num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            WriteASM_S(sw, shares, TT);
            WriteASM_P(sw, shares);
            sw.Close();
            fs.Close();
        }
        //用ASM写出软件实现TI的代码
        //用sharesgroup分组
        public void Print_TI1b_ASM_SharesGroup(string path, long num, string sname, int shares, byte[] TT)
        {

            //Write the ASM
            String filename = String.Format(path + sname + "{0}_R{1}_{2}.S", size, round, num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            WriteASM_S_SharesGroup(sw, shares, TT);
            //WriteASM_P(sw, shares);
            sw.Close();
            fs.Close();
        }
        //用C写出软件实现TI的代码
        public void Print_TI1b_C(string path, long num, string sname, int shares, byte[] TT)
        {
            int[] ANF = MoebiusTrans(size * (shares - 1), TT);
            //Compute ANF terms
            ANFterm[] anf = ConvertANF(ANF, shares);
            //Write the ASM
            String filename = String.Format(path + sname + "{0}_R{1}_{2}.C", size, round, num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            int len_r = size * (shares - 1);
            int len_e = size * (shares + 1);
            //r0 is the input register, contains shared data like x00x01x02x00x10x11x12x10....
            //Compute the shift version
            //Use r{len_r} as temp
            //Use r{len_r+2} as mask:0x0000ffff
            //Use r{len_r+1} as mask:0b'0001000100010001
            //Use r{len_r+3} as mask:0b'0111011101110111
            sw.WriteLine("   unsigned int reg[{0}];",len_r+5);
            sw.WriteLine("   reg[0] = input;");
            sw.WriteLine("   reg[{0}] = GetMask;",len_r+1);
            sw.WriteLine("   reg[{0}] = ZeroMask;", len_r + 3); 
            sw.WriteLine("   reg[{0}]= AndMask;", len_r + 2);
            sw.WriteLine("   reg[{0}]= PMask;", len_r + 4);
            for (int i = 1; i <len_r; i++)
            {
                if ((i % (shares - 1) == 0) && (i > 1))
                {
                    //Right shift varables=Rotated shift by 4
                    sw.WriteLine("   reg[{0}]=reg[{1}]>>4;", i, i - 2);//reg[i]=reg[i-2]>>4;
                    sw.WriteLine("   reg[{0}]=reg[{1}]<<{2};", len_r, i - 2, len_e - 4);//temp=reg[i-2]<<12
                    sw.WriteLine("   reg[{0}]=reg[{1}]^reg[{0}];", i, len_r);//reg[i]=reg[i]^temp;
                    sw.WriteLine("   reg[{0}]=reg[{0}]&reg[{1}];", i, len_r + 2);//reg[i]=reg[i]&mask;
                }
                else
                {
                    //Right shift shares by 1=Rotated shift shares by 1
                    sw.WriteLine("   reg[{0}]=reg[{1}]>>1;", i, i - 1);//reg[i]=reg[i-1]>>1
                    sw.WriteLine("   reg[{0}]=reg[{1}]&reg[{2}];", len_r, i , len_r+1);//temp=reg[i]&GetMask
                    sw.WriteLine("   reg[{0}]=reg[{1}]&reg[{2}];", i, i, len_r + 3);//reg[i]=reg[i]&ZeroMask
                    sw.WriteLine("   reg[{0}]=reg[{1}]<<3;", len_r, len_r);//temp=temp<<3
                    sw.WriteLine("   reg[{0}]=reg[{1}]^reg[{0}];", i, len_r);//reg[i]=reg[i]^temp
                }
            }
            //Sort anf as degree and delta_index
            Array.Sort(anf);
            //Use r{len_r} as result
            //Use r{len_r+1} as temp
            sw.WriteLine("   reg[{0}]=0;", len_r);
            for (int i = 0; i < anf.Length; i++)
            {
                if (anf[i].d == 1)//linear term
                {
                    sw.WriteLine("   reg[{0}]=reg[{0}]^reg[{1}];", len_r, anf[i].base_index);
                }
                else//only deal with degree 2 funcitons
                {
                    sw.WriteLine("   reg[{0}]=reg[{1}]&reg[{2}];", len_r + 1, anf[i].base_index, anf[i].base_index+anf[i].other_index[0]);//temp1=r0&r_delta
                    //Xor temp to result
                    sw.WriteLine("   reg[{0}]=reg[{0}]^reg[{1}];", len_r, len_r + 1);

                }
            }
            //Do Pmatrix
            sw.WriteLine("   reg[{1}]=reg[{0}]<<4;", len_r, len_r + 1);//temp1=result<<4;
            sw.WriteLine("   reg[{0}]=reg[{0}]&AndMask;", len_r + 1);//temp1=temp1&AndMask;
            sw.WriteLine("   reg[{0}]=reg[{0}]>>{1};", len_r, len_e-4);//result=result>>12;
            sw.WriteLine("   reg[{0}]=reg[{0}]&0x07;", len_r);//result=result&0x07;
            sw.WriteLine("   reg[{0}]=PMTable[reg[{0}]];",len_r);//result=PMTable[resullt];
            sw.WriteLine("   reg[{0}]=reg[{1}]^reg[{0}];",len_r, len_r + 1);//result=temp1^result;
            //Do Pmatrix2
            sw.WriteLine("   reg[{1}]=reg[{0}]<<4;", len_r, len_r + 1);//temp1=result<<4;
            sw.WriteLine("   reg[{0}]=reg[{0}]&AndMask;", len_r + 1);//temp1=temp1&AndMask;
            sw.WriteLine("   reg[{0}]=reg[{0}]>>{1};", len_r, len_e - 4);//result=result>>12;
            sw.WriteLine("   reg[{0}]=reg[{0}]&0x07;", len_r);//result=result&0x07;
            sw.WriteLine("   reg[{0}]=PMTable[reg[{0}]];", len_r);//result=PMTable[resullt];
            sw.WriteLine("   reg[{0}]=reg[{1}]^reg[{0}];", len_r, len_r + 1);//result=temp1^result;

            sw.WriteLine("   return reg[{0}];", len_r);

            sw.Close();
            fs.Close();
        }
        //用C写出软件实现TI的代码
        //输入按shares分组
        public void Print_TI1b_C_SharesGroup(string path, long num, string sname, int shares, byte[] TT)
        {
            int[] ANF = MoebiusTrans(size * (shares - 1), TT);
            //Compute ANF terms
            ANFterm[] anf = ConvertANF(ANF, shares);
            //Write the ASM
            String filename = String.Format(path + sname + "{0}_R{1}_{2}.C", size, round, num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            int len_r = size * (shares - 1);
            int len_e = 2*size * shares;
            //r0 is the input register, contains shared data like x00x10x20x30x01x11x21x31....
            //Compute the shift version
            //Use r{len_r} as temp
            //Use r{len_r+2} as mask:0x00ffffff
            //Use r{len_r+1} as mask:0x000f0f0f
            sw.WriteLine("   unsigned int reg[{0}];", len_r + 5);
            sw.WriteLine("   reg[0] = input;");
            sw.WriteLine("   reg[{0}] = GetMask;", len_r + 1);
            sw.WriteLine("   reg[{0}]= AndMask;", len_r + 2);
            sw.WriteLine("   reg[{0}]= PMask;", len_r + 4);
            for (int i = 1; i < len_r; i++)
            {
                if (i % (size) == 0)//shift shares
                {
                    //Right shift shares==shift 2*size
                    sw.WriteLine("   reg[{0}]=reg[{1}]>>{2};", i, i - size,2*size);//reg[i]=reg[i-size]>>8;
                    sw.WriteLine("   reg[{0}]=reg[{1}]<<{2};", len_r, i - size, len_e - 2 * size);//temp=reg[i-size]<<16
                    sw.WriteLine("   reg[{0}]=reg[{1}]^reg[{0}];", i, len_r);//reg[i]=reg[i]^temp;
                    sw.WriteLine("   reg[{0}]=reg[{0}]&reg[{1}];", i, len_r + 2);//reg[i]=reg[i]&Andmask;
                }
                else//Shift bits
                {
                    //Right shift bits by 1
                    sw.WriteLine("   reg[{0}]=reg[{1}]>>1;", i, i - 1);//reg[i]=reg[i-1]>>1
                    sw.WriteLine("   reg[{0}]=reg[{1}]&reg[{2}];", len_r, i, len_r + 1);//temp=reg[i]&GetMask
                    sw.WriteLine("   reg[{0}]=reg[{1}]<<{2};", len_r, len_r,size);//temp=temp<<4
                    sw.WriteLine("   reg[{0}]=reg[{1}]&reg[{2}];", i, i, len_r + 1);//reg[i]=reg[i]&GetMask
                    sw.WriteLine("   reg[{0}]=reg[{1}]^reg[{0}];", i, len_r);//reg[i]=reg[i]^temp
                }
            }
            //Sort anf as degree and delta_index
            Array.Sort(anf);
            //Use r{len_r} as result
            //Use r{len_r+1} as temp
            sw.WriteLine("   reg[{0}]=0;", len_r);
            for (int i = 0; i < anf.Length; i++)
            {
                if (anf[i].d == 1)//linear term
                {
                    sw.WriteLine("   reg[{0}]=reg[{0}]^reg[{1}];", len_r, anf[i].base_index);
                }
                else//only deal with degree 2 funcitons
                {
                    sw.WriteLine("   reg[{0}]=reg[{1}]&reg[{2}];", len_r + 1, anf[i].base_index, anf[i].base_index + anf[i].other_index[0]);//temp1=r0&r_delta
                    //Xor temp to result
                    sw.WriteLine("   reg[{0}]=reg[{0}]^reg[{1}];", len_r, len_r + 1);

                }
            }
         

            sw.WriteLine("   return reg[{0}];", len_r);

            sw.Close();
            fs.Close();
        }
        public void PrintScript(StreamWriter swScript,string pathname,string Sname,long num)
        {
            swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/"+pathname+"/"+String.Format(Sname+"{0}_R{1}_{2}.v",size,round,num)+"\"}");
            swScript.WriteLine("compile -exact_map");
            swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/"+pathname+"/areareport_" + num + ".txt\"");
            swScript.WriteLine("remove_design -designs");
        }
        //验证TITable还原后能否得到原始表
        public bool VerifyTITable(int[] TITT, int[] Stable, int shares)
        {
            int[] NewS = new int[Stable.Length];
            int full = shares * size;
            if (!CheckPermutaion(Stable))
            {
                System.Console.WriteLine("Sbox not a permutation!");
                return false;
            };
            bool[] Used = new bool[(0x1 << full)];
            for (int i = 0; i < (0x1 << full); i++)
            {
                int[] sharedinput = new int[size];
                int[] sharedoutput = new int[size];
                //Get input as an array
                int mask = (0x1 << (shares))-1;
                int temp=i;
                for (int j = 0; j < size; j++)
                {
                    sharedinput[size - 1 - j] = temp & mask;
                    temp = temp >> shares;
                }
                //Compute for each share
                for (int j = 0; j < shares; j++)
                {
                    for(int k=0;k<size;k++)
                    {
                         //Get input
                           int t = 0;
                            int newmask=(0x1<<(shares-1))-1;
                            for (int m = 0; m < size; m++)
                            {
                                t = t | (sharedinput[(m+size-k)%size] & newmask);
                                //t = t | (sharedinput[(m+k)%size] & newmask);
                                if(m!=size-1)
                                    t = t << (shares - 1);
                            }
                          //Compute output
                            sharedoutput[k] ^= TITT[t];
                            if (j != shares - 1)
                                sharedoutput[k] = sharedoutput[k] << 1;
                    }
                    //Shift Shares
                    for (int k = 0; k < size; k++)
                    {
                        sharedinput[k] = ((sharedinput[k]&0x01)<<(shares-1)) | (sharedinput[k] >> 1);
                    }

                }
                //Mark as Used
                int tempv = 0;
                for (int j = 0; j < size; j++)
                {
                    tempv = tempv ^ sharedoutput[j];
                    if (j != size - 1)
                        tempv = tempv << shares;
                }
                if (Used[tempv] == true)
                {
                    System.Console.WriteLine("Not an TI Permutation!");
                    return false;
                }
                Used[tempv] =true;
                //Get Result
                int input = 0;
                int output = 0;
                for (int j = 0; j < size; j++)
                {
                    if (HW(sharedinput[j], shares) % 2 == 1)
                        sharedinput[j] = 1;
                    else
                        sharedinput[j] = 0;
                    if (HW(sharedoutput[j], shares) % 2 == 1)
                        sharedoutput[j] = 1;
                    else
                        sharedoutput[j] = 0;
                    input = (input ^ sharedinput[j]);
                    //output = (output ^ sharedoutput[j]);
                    if (j != size - 1)
                    {
                        input = input << 1;
                    //    output = output << 1;
                    }
                    if (sharedoutput[j] == 1)
                        output = SetBit(output, j);
                    //if (sharedinput[j] == 1)
                    //    input = SetBit(input, j);
                }
                if (Stable[input] != output)
                    return false;
            }
            
            return true;
        }
        //验证TITable还原后能否得到原始表
        //按同一个shares放在一起的方式处理
        public bool VerifyTITable_SharesGroup(int[] TITT, int[] Stable, int shares)
        {
            int[] NewS = new int[Stable.Length];
            int full = shares * size;
            if (!CheckPermutaion(Stable))
            {
                System.Console.WriteLine("Sbox not a permutation!");
                return false;
            };
            bool[] Used = new bool[(0x1 << full)];
            for (int i = 0; i < (0x1 << full); i++)
            {
                int[] sharedinput = new int[shares];
                int[] sharedoutput = new int[shares];
                //Get input as an array
                int mask = (0x1 << (size)) - 1;
                int temp = i;
                for (int j = 0; j < shares; j++)
                {
                    //sharedinput[shares-1-j] = temp & mask;
                    sharedinput[shares - 1 - j] = temp & mask;
                    temp = temp >> size;
                }
                //Compute for each share
                for (int j = 0; j < shares; j++)
                {
                    sharedoutput[j] = 0;
                    //Get each bit
                    for (int k = 0; k < size; k++)
                    {
                        //Get input
                        int t = 0;
                        for (int m = 0; m < shares-1; m++)
                        {
                            t = t ^ sharedinput[(m + j+shares-2) % shares];
                            if (m != shares - 2)
                                t = t << size;
                        }
                        //Compute output
                        if (TITT[t] == 1)
                            sharedoutput[j] = SetBit(sharedoutput[j], k);
                        //sharedoutput[j] ^= TITT[t];
                        //if (k != size - 1)
                        //{
                        //    sharedoutput[j] = sharedoutput[j] << 1;
                        //}
                        //Shift bits
                        //right shift bits by 1
                        for (int m = 0; m < shares; m++)
                        {
                            sharedinput[m] = (sharedinput[m] >> 1) ^ ((sharedinput[m] & 0x01) << (size - 1));
                        }
                    }
                }
                //Mark as Used
                int tempv = 0;
                for (int j = 0; j < shares; j++)
                {
                    tempv = tempv ^ sharedoutput[j];
                    if (j != shares - 1)
                        tempv = tempv << size;
                }
                if (Used[tempv] == true)
                {
                    System.Console.WriteLine("Not an TI Permutation!");
                    return false;
                }
                Used[tempv] = true;
                //Get Result
                int input = 0;
                int output = 0;
                for (int j = 0; j < shares; j++)
                {
                   
                    input = (input ^ sharedinput[j]);
                    output = (output ^ sharedoutput[j]);
                }
                if (Stable[input] != output)
                {
                    System.Console.WriteLine("Error!");
                   // return false;
                }
            }

            return true;
        }
        #endregion



        //SITI+C
        public void OneRoundTrans_SITIC(int[] table, int[] S, int C)
        {

            for (int i = 0; i < len; i++)
            {
                table[i] = S[table[i]] ^ C;
            }
        }

        //SITI+*2P(左移位)
        public void OneRoundTrans_SITIM2P(int[] table, int[] S, int Pno)
        {

            for (int i = 0; i < len; i++)
            {
                table[i] = S[table[i]];
                if ((table[i] >> (size - 1)) > 0)
                    table[i] = ((table[i] << 1)&((0x1<<size)-1)) ^ ((Pno << 1) ^ 0x01);
                else
                    table[i] = ((table[i] << 1) & ((0x1 << size) - 1));
            }
        }
        //SITI+*4P
        public void OneRoundTrans_SITIM4P(int[] table, int[] S, int Pno)
        {

            for (int i = 0; i < len; i++)
            {
                table[i] = S[table[i]];
                if ((table[i] >> (size - 1)) > 0)
                    table[i] = ((table[i] << 1) & ((0x1 << size) - 1)) ^ ((Pno << 1) ^ 0x01);
                else
                    table[i] = ((table[i] << 1) & ((0x1 << size) - 1));
                if ((table[i] >> (size - 1)) > 0)
                    table[i] = ((table[i] << 1) & ((0x1 << size) - 1)) ^ ((Pno << 1) ^ 0x01);
                else
                    table[i] = ((table[i] << 1) & ((0x1 << size) - 1));
            }
        }
        //SITI+MatrixP
        public void OneRoundTrans_SITIMatrixP(int[] table, int[] S, int FieldNo,int MultiNo)
        {

            for (int i = 0; i < len; i++)
            {
                table[i] = S[table[i]];
                table[i] = FieldMultiplication(table[i], FieldNo, MultiNo);
            }
        }
        //SITI+*2P(右移位)
        public void OneRoundTrans_SITIM2P_Right(int[] table, int[] S, int Pno)
        {

            for (int i = 0; i < len; i++)
            {
                table[i] = S[table[i]];
                if ((table[i]&0x01) > 0)
                    table[i] = ((table[i] >> 1) ) ^ ((Pno) ^ (0x1<<(size-1)));
                else
                    table[i] = ((table[i] >> 1) );
            }
        }
        //SITI+AnyP
        public void OneRoundTrans_SITIAnyP(int[] table, int[] S, int[][] PMatrix)
        {

            for (int i = 0; i < len; i++)
            {
                table[i] = S[table[i]];
                table[i] = LinearP(PMatrix, table[i], size);
            }
        }
        //SITI 4 bit Sbox搜索,添加常数
        //Results: 存在15-15的定差分，无论增加多少轮，如何增加常数或者增加线性变换，都不能破坏定差分
        public void SearchOptimal_SIConstant_4(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的4bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, 4);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            long length = Stable.Length*16 ;
            for (long num =0; num < length; num++)
            {
                if ((num&0xff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0*num / length);
                }
                //从num中获取S的序号
                int Sno = (int)(num /16);
                int C = (int)(num % 16);
                //int C = 1;
                //int Sno = (int)num;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIC(table,Stable[Sno],C);
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

        }


        //SITI 4 bit Sbox搜索,添加线性层(*2线性层)
        //Results: 4轮能找到16个最优的
        public void SearchOptimal_SIP_4(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的4bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, 4);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            int Psize=0x1<<(size-1);
            long length = Stable.Length * (Psize);
            for (long num = 0; num < length; num++)
            {
                if ((num & 0xff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0 * num / length);
                }
                //从num中获取S的序号
                int Sno = (int)(num / Psize);
                int Pno = (int)(num % Psize);
                int[][] PMatrix=BuildM2PMatrix(Pno, size);
                //检验可逆性
                if (ReverseM(PMatrix, size) == 0)//不可逆
                    continue;
                //int C = 1;
                //int Sno = (int)num;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIM2P_Right(table, Stable[Sno],Pno);
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
                        for (int j = 0; j <(0x1<<size); j++)
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

        }
        //SITI 5 bit Sbox搜索,添加常数
        //Results: 一部分为32的定差分，另一半有16的定差分
        public void SearchOptimal_SIConstant_5(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的5bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, 5);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            long length = Stable.Length*32;
            for (long num = 0; num < length; num++)
            {
                if ((num & 0xff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0 * num / length);
                }
                //从num中获取S的序号
                int Sno = (int)(num /32);
                int C = (int)(num % 32);
                //int C = 1;
                //int Sno = (int)num;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIC(table, Stable[Sno], C);
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

        }

        //SITI 8 bit Sbox搜索
        //Results: 
        public void SearchOptimal_SI_8(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的5bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, 8);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            long length = Stable.Length;
            for (long num = 0; num < length; num++)
            {
                if ((num & 0xff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0 * num / length);
                }
                //从num中获取S的序号
                int Sno = (int)(num);
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIC(table, Stable[Sno], 0);
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

        }


        //SITI 5 bit Sbox搜索,添加线性层(*2线性层)
        //Results: 4轮差分为6有144个
        public void SearchOptimal_SIP_5(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的4bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, size);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            int Psize = 0x1 << (size - 1);
            long length = Stable.Length * (Psize);
            for (long num = 0; num < length; num++)
            {
                if ((num & 0xff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0 * num / length);
                }
                //从num中获取S的序号
                int Sno = (int)(num / Psize);
                int Pno = (int)(num % Psize);
                int[][] PMatrix = BuildM2PMatrix(Pno, size);
                //检验可逆性
                if (ReverseM(PMatrix, size) == 0)//不可逆
                    continue;
                //int C = 1;
                //int Sno = (int)num;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIM2P(table, Stable[Sno], Pno);
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
                        for (int j = 0; j < (0x1 << size); j++)
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

        }
        //SITI 任意bit Sbox搜索,添加线性层(*2线性层)
        public void SearchOptimal_SIP(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的4bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, size);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            int Psize = 0x1 << (size - 1);
            long length = Stable.Length * (Psize);
            for (long num = 0; num < length; num++)
            {
                if ((num & 0xff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0 * num / length);
                }
                //从num中获取S的序号
                int Sno = (int)(num / Psize);
                int Pno = (int)(num % Psize);
                int[][] PMatrix = BuildM2PMatrix(Pno, size);
                //检验可逆性
                if (ReverseM(PMatrix, size) == 0)//不可逆
                    continue;
                //int C = 1;
                //int Sno = (int)num;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIM2P(table, Stable[Sno], Pno);
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
                        for (int j = 0; j < (0x1 << size); j++)
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

        }

        //SITI 任意bit Sbox搜索,添加线性层(*4线性层)
        public void SearchOptimal_SIM4P(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear,string path, string scriptfilename,string pathname,string sname)
        {
            int shares = 3;
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的4bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, size);
            byte[][] TT = ReadShiftInvariantTT_OneBit(Sbinfile, size);

            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);

            FileStream fsScript = new FileStream(path+scriptfilename, FileMode.Create);
            StreamWriter swScript = new StreamWriter(fsScript);

            //写Script头
            swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            swScript.WriteLine("set synthetic_library generic.sdb");
            swScript.WriteLine("set target_library \"fast.db slow.db\"");
            swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");

            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            int Psize = 0x1 << (size - 1);
            long length = Stable.Length * (Psize);
             Parallel.For(0, length, new ParallelOptions { MaxDegreeOfParallelism =8 }, num =>
            {
                if ((num & 0xffff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%,mindiff={2}", num, 100.0 * num / length,mindiff);
                }
                //从num中获取S的序号
                int Sno = (int)(num / Psize);
                int Pno = (int)(num % Psize);
                int[][] PMatrix = BuildM2PMatrix(Pno, size);
                //检验可逆性
                if (ReverseM(PMatrix, size) == 0)//不可逆
                    return;
                //int C = 1;
                //int Sno = (int)num;
                //table初始化为恒等变换
                int[] table = new int[len];
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIM4P(table, Stable[Sno], Pno);
                }
                if (!CheckPermutaion(table))
                {
                    System.Console.WriteLine("Error!");
                    return;
                }
                    
                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, size);

                if (diff < mindiff)
                    Interlocked.Exchange(ref mindiff,diff);

                if ((num & 0xffff) == 0)
                    System.Console.WriteLine("diff={0},num={1:x}", diff, num);

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
                        
                        System.Console.WriteLine("*********diff={0},nL={1},num={2}", diff, nonLinear, num);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},Sno={2:x},Pno={3:x}", diff, nonLinear, Sno,Pno);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < (0x1 << size); j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("内部函数={0}", num);

                        //Get TI bit
                        int[] ANF = MoebiusTrans(size, TT[Sno]);//Get ANF
                        //Compute the TI form from ANF
                        TI_wrapper tw = new TI_wrapper(size, shares, ANF);
                        tw.Compute_From_ANF();
                        int[] TITT = tw.Get_TruthTable_F1();
                        if (!VerifyTITable(TITT, Stable[Sno], shares))
                        {
                            System.Console.WriteLine("TI Sbox results Incorrect!");
                            return;
                        }


                        byte[] TITTb=TTTransformation(TITT);
                         //Write ASM TI implementation
                       // Print_TI1b_ASM(path, num, sname, 3, TITTb);
                        //Print_TI1b_ASM(path, num, sname, 3, TITTb);
                       
                        //评估代价
                        WriteSIM4PCost_TI1b(sw, TITTb, Pno);
                        WriteSIM4PCost_Unprotected(sw, TT[Sno], Pno);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();
                        //写出Sbox
                        //Print_FullSbox_Unprotected(path, table, num,sname);
                        //Print_SI1bit_Unprotected(path, Stable[Sno], num, sname);
                         Print_TI1b(path, TITT, num,sname,3);
                        //写出脚本
                        PrintScript(swScript,pathname, sname, num);
                    }

                    }

                }
            });
            System.Console.WriteLine("number of possible solutions: {0}", count);
            System.Console.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.WriteLine("mindiff={0},max non-linear={1}", mindiff, maxnL);
            sw.Close();
            fs.Close();
            swScript.Close();
            fsScript.Close();
        }
        //SITI 任意bit Sbox搜索,添加任意线性层
        public void SearchOptimal_SIAnyP(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的4bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, size);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            long Psize = 0x1 << (size*size);
            long length = Stable.Length * (Psize);
            for (long num = 0; num < length; num++)
            {
                if ((num & 0xffff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0 * num / length);
                }
                //从num中获取S的序号
                long Pno = (num / Stable.Length);
                int[][] PMatrix = BuildAnyPMatrix(Pno, size);
                //检验可逆性
                if (ReverseM(PMatrix, size) == 0)//不可逆
                    continue;

                int Sno = (int)(num % Stable.Length);
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIAnyP(table, Stable[Sno], PMatrix);
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
                        for (int j = 0; j < (0x1 << size); j++)
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

        }
        //SITI 任意bit Sbox搜索,添加线性层(有限域乘法)
        public void SearchOptimal_SIMatrixP(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的4bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, size);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            long Psize = (0x1 << (size - 1))*((0x1<<size)-2);//前一个是选择P的域，后一个选择域元素，去除0和1
            long length = Stable.Length * (Psize);
            for (long num = 0; num < length; num++)
            {
                if ((num & 0xffff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%", num, 100.0 * num / length);
                }
                //从num中获取S的序号
                int Sno = (int)(num % Stable.Length);
                long Pno = (num / Stable.Length);
                int FieldPno =(int)( Pno % (0x1 << (size - 1)));//域选择
                int MultiPno = (int) (Pno / (0x1 << (size - 1)));//域元素
                MultiPno = MultiPno + 2;
                int[][] PMatrix = BuildM2PMatrix(FieldPno, size);
                //检验可逆性
                if (ReverseM(PMatrix, size) == 0)//不可逆
                    continue;
                //int C = 1;
                //int Sno = (int)num;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIMatrixP(table, Stable[Sno], FieldPno, MultiPno);
                }
                if (!CheckPermutaion(table))
                    continue;
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
                        System.Console.WriteLine("*********diff={0},nL={1},num={2:x},multi={3:x}", diff, nonLinear, num,MultiPno);
                        count++;
                        sw.WriteLine("******************");
                        sw.WriteLine("diff={0},nL={1},num={2:x}", diff, nonLinear, num);
                        sw.WriteLine("输出表：");
                        for (int j = 0; j < (0x1 << size); j++)
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

        }

        //SITI 8 bit Sbox搜索,添加线性层(*2线性层)
        //Results: 
        public void SearchOptimal_SIP_8(string Sbinfile, string filename, int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            long[] item = new long[round];
            int mindiff = 256;
            int maxnL = 0;
            //读取所有有直接TI的4bit SI 2次置换
            int[][] Stable = ReadShiftInvariantTT_Sbox(Sbinfile, size);


            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine("round={0},size={1},shift={2}", round, size, shift);
            int Psize = 0x1 << (size - 1);
            long length = Stable.Length * (Psize);
            for (long num = 0; num < length; num++)
            {
                if ((num & 0xffff) == 0)
                {
                    System.Console.WriteLine("num={0:x},percent={1}%,mindiff={2}", num, 100.0 * num / length,mindiff);
                }
                //从num中获取S的序号
                int Sno = (int)(num / Psize);
                int Pno = (int)(num % Psize);
                int[][] PMatrix = BuildM2PMatrix(Pno, size);
                //检验可逆性
                if (ReverseM(PMatrix, size) == 0)//不可逆
                    continue;
                //int C = 1;
                //int Sno = (int)num;
                //table初始化为恒等变换
                for (int i = 0; i < len; i++)
                {
                    table[i] = i;
                }

                for (int i = 0; i < round; i++)
                {
                    OneRoundTrans_SITIM2P(table, Stable[Sno], Pno);
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
                        for (int j = 0; j < (0x1 << size); j++)
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

        }
    }


}
