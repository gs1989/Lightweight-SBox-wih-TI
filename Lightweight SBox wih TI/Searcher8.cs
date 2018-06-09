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
        //cost[0]:AND的数量;cost[1]:XOR的数量;cost[2]:Shift的数量;cost[3]:M0 Cycle总和
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
                        {
                            Cost[2] = Cost[2] + 2;
                            Cost[3] = Cost[3] + 8;//Mov+Shift+Mov+Shift+AND+EOR
                        }
                        else//有一个是x0
                        {
                            Cost[2] = Cost[2] + 1;
                            Cost[3] = Cost[3] + 6;//Mov+Mov+Shift+AND+EOR
                        }
                    }
                    if (HW(i, size) == 1)//1次项
                    {
                        if ((i & 0x01) == 0)//不是x0
                        {
                            Cost[2] = Cost[2] + 1;
                            Cost[3] = Cost[3] + 4;//Mov+Shift+EOR
                        }
                        else
                            Cost[3] = Cost[3] + 2;//Mov+EOR

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
            
            sw.WriteLine("Unprotected Implementation:");
            sw.WriteLine("Software Cost: ARM M0");
            int[] SoftwareCost = ANFCost_Software(TT, size);
            sw.WriteLine("Software Cycles of S: AND={0},XOR={1},RotatedShift={2},Cycles={3}", SoftwareCost[0] * round, SoftwareCost[1] * round, SoftwareCost[2] * round, SoftwareCost[3] * round);
            int ScostM0 = SoftwareCost[3] * round;
            int ScostM3 = (SoftwareCost[3] - SoftwareCost[2] * 2) * round;
            int sliceblock = 32 / size;
            int Pcost = (6 * sliceblock + 10) * (round - 1);
            int AllCostM0 = (ScostM0 + 2 * Pcost);
            int AllCostM3 = (ScostM3 + 2 * Pcost);
            sw.WriteLine("Software Cycles of P: Cycles={0}", 2 * Pcost);
            sw.WriteLine("Overall cycles={0}, Average={1}", AllCostM0, AllCostM0 / (double)sliceblock);

            sw.WriteLine("Software Cost: ARM M3");
            sw.WriteLine("Software Cycles of S: AND={0},XOR={1},RotatedShift={2},Cycles={3}", SoftwareCost[0] * round, SoftwareCost[1] * round, SoftwareCost[2] * round, ScostM3);
            sw.WriteLine("Software Cycles of P: Cycles={0}", 2 * Pcost);
            sw.WriteLine("Overall cycles={0}, Average={1}", AllCostM3, AllCostM3 / (double)sliceblock);

            sw.Flush();
        }
        //评估整体unprotecked Sbox的实现代价
        public void WriteSIM4PSoftwareCost_TI1b(StreamWriter sw, byte[] TT1b, int Pno)
        {
            int shares = 3;
            sw.WriteLine("Software Cost: ARM M0");

            int[] SoftwareCost = ANFCost_Software(TT1b, size * (shares - 1));
            sw.WriteLine("Software Cycles of S: AND={0},XOR={1},RotatedShift={2},Cycles={3}", SoftwareCost[0] * round, SoftwareCost[1] * round, SoftwareCost[2] * round, SoftwareCost[3] * round);
            int ScostM0 = SoftwareCost[3]*round;
            int ScostM3 = (SoftwareCost[3] - SoftwareCost[2] * 2) * round;
            int sliceblock = 32 / size;
            int Pcost = (6 * sliceblock + 10) * (round - 1);
            int AllCostM0 = shares * (ScostM0 + 2 * Pcost);
            int AllCostM3 = shares * (ScostM3 + 2 * Pcost);
            sw.WriteLine("Software Cycles of P: Cycles={0}", 2 * Pcost);
            sw.WriteLine("Overall cycles={0}, Average={1}", AllCostM0, AllCostM0 / (double)sliceblock);
            
            sw.WriteLine("Software Cost: ARM M3");
            sw.WriteLine("Software Cycles of S: AND={0},XOR={1},RotatedShift={2},Cycles={3}", SoftwareCost[0] * round, SoftwareCost[1] * round, SoftwareCost[2] * round, ScostM3);
            sw.WriteLine("Software Cycles of P: Cycles={0}", 2 * Pcost);
            sw.WriteLine("Overall cycles={0}, Average={1}", AllCostM3, AllCostM3 / (double)sliceblock);

            sw.Flush();
        }
        //用txt写出S表
        public void Print_FullSbox_Unprotected(string path, int[] Sbox,long num,string sname)
        {
            String filename=String.Format( path+sname+"{0}_R{1}_{2}.txt",size,round,num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            

            for (int i = 0; i < Sbox.Length; i++)
                if(i!=Sbox.Length-1)
                   sw.Write("{0},", Sbox[i]);
                else
                    sw.Write("{0}", Sbox[i]);

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
        //2018.5.10 Shift all shares at once
        public void WriteASM_S_8(StreamWriter sw, int shares, byte[] TT)
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
            //Pass r1 as memory space that contains x>>>0 to x>>>15
            //Pass r2 as ShareMask2=0xcccc...c
            //Pass r3 as ShareMask1=0x2222...2
            //Use r4 as temp1
            //Use r5 as constant 4
            //Use r6 as AndMask (if needed)
            //Use r7 as temp
            if (len_e < 32)//Create AndMask
            {
                sw.WriteLine("   movs r6,#1");
                sw.WriteLine("   lsls r6,#{0}",len_e);
                sw.WriteLine("   subs r6,#1", len_e);
            }
                //Store x0
                sw.WriteLine("   str r0,[r1,#0]");
                //Compute Shift 1
                //Higher two shares
                sw.WriteLine("   mov r7,r0");//r7=x
                sw.WriteLine("   ands r7,r2");//r7=x&ShareMask2
                sw.WriteLine("   lsrs r7,#1");//temp=r0>>1;
                //Lower one share
                sw.WriteLine("   mov r4,r0");//r4=x
                sw.WriteLine("   ands r4,r3");//r4=x&ShareMask1
                sw.WriteLine("   lsrs r4,#1");//r4=(x&ShareMask1)>>1
                sw.WriteLine("   lsls r5,r4,#3");//r5=(x&ShareMask1)<<2
                sw.WriteLine("   eors r5,r4");//r5=((x&ShareMask1)<<2)^((x&ShareMask1)>>1)
                //Paste together
                sw.WriteLine("   eors r7,r5");//r7=(x>>>1share)
               //Store x1
                sw.WriteLine("   nop");
                sw.WriteLine("   str r7,[r1,#4]");
               //Set constant
                sw.WriteLine("   movs r5,#4");

                //Start from i=2
                for (int i = 2; i < len_r; i++)
                {
                    if ((i % (shares - 1) == 0))
                    {
                        //Shift r0 by 4
                        if (len_e == 32)
                            sw.WriteLine("   rors r0,r5");//r0=r0>>>4
                        else//Do not have ror
                        {
                            sw.WriteLine("   lsls r4,r0,#12");//r4=r0<<12
                            sw.WriteLine("   ands r4,r6");//r4=(r0<<12)&AndMask
                            sw.WriteLine("   lsrs r0,r5");//r0=r0>>4
                            sw.WriteLine("   eors r0,r4");//r0=((r0<<12)&AndMask)^(r0>>4)=(r0>>>4)
                        }
                        //Store x[i*4]
                        sw.WriteLine("   nop");
                        sw.WriteLine("   str r0,[r1,#{0}]", i * 4);
                    }
                    else
                    {
                        //Shift r7 by 4
                        if (len_e == 32)
                            sw.WriteLine("   rors r7,r5");//r7=r7>>>4
                        else//Do not have ror
                        {
                            sw.WriteLine("   movs r4,#0");//Clear temp1
                            sw.WriteLine("   lsls r4,r7,#12");//r4=r7<<12
                            sw.WriteLine("   ands r4,r6");//r4=(r7<<12)&AndMask
                            sw.WriteLine("   lsrs r7,r5");//r7=r7>>4
                            sw.WriteLine("   eors r7,r4");//r7=((r7<<12)&AndMask)^(r7>>4)=(r7>>>4)
                        }
                        //Store x[i*4]
                        sw.WriteLine("   nop");
                        sw.WriteLine("   str r7,[r1,#{0}]", i * 4);
                    }
                }
            

            //Sort anf as degree and delta_index
            Array.Sort(anf);
            //Use r0 as result
            //Use r2 as temp1
            //Use r3 as temp2
            sw.WriteLine("   movs r0,#0");
            for (int i = 0; i < anf.Length; i++)
            {

                if (anf[i].d == 1)//linear term
                {
                    //Clear Op2 before another LDR
                    sw.WriteLine("   nop");
                    sw.WriteLine("   ldr r2,[r1,#{0}]", anf[i].base_index * 4);//temp1=reg[];
                   
                   
                    //if ((anf[i].base_index & 0x01) == 0)
                    //    sw.WriteLine("   mov  r2,r4");//r2=x
                    //else
                    //    sw.WriteLine("   mov  r2,r5");//r2=x>>>1
                    //sw.WriteLine("   movs r6,#{0}", 4*anf[i].base_index / 2);
                    //sw.WriteLine("   rors r2,r6");//r2>>>anf[i].base_index
                }
                else//only deal with degree 2 funcitons
                {
                    //Clear Op2 before another LDR
                    sw.WriteLine("   nop");
                    sw.WriteLine("   ldr r2,[r1,#{0}]", anf[i].base_index * 4);//temp1=reg[];
                    //Clear Op2 before another LDR
                    sw.WriteLine("   nop");
                    sw.WriteLine("   ldr r3,[r1,#{0}]", (anf[i].base_index + anf[i].other_index[0]) * 4);//temp2=reg[];
                    sw.WriteLine("   ands r2,r3");//temp1=temp1&temp2
                    //if ((anf[i].base_index & 0x01) == 0)
                    //    sw.WriteLine("   mov  r2,r4");//r2=x
                    //else
                    //    sw.WriteLine("   mov  r2,r5");//r2=x>>>1
                    //sw.WriteLine("   movs r6,#{0}", 4 * (anf[i].base_index / 2));
                    //sw.WriteLine("   rors r2,r6");//r2>>>anf[i].base_index
                    //if (((anf[i].base_index + anf[i].other_index[0]) & 0x01) == 0)
                    //    sw.WriteLine("   mov  r3,r4");//r3=x
                    //else
                    //    sw.WriteLine("   mov  r3,r5");//r3=x>>>1
                    //sw.WriteLine("   movs r6,#{0}", 4 * ((anf[i].base_index + anf[i].other_index[0]) / 2));
                    //sw.WriteLine("   rors r3,r6");//r3>>>(anf[i].base_index + anf[i].other_index[0])
                    ////And
                    //sw.WriteLine("   ands r2,r3");//temp1=temp1&temp2
                }
                //Xor temp to result
                sw.WriteLine("   eors r0,r2");
            }
            for (int i = 0; i < 8; i++)
                sw.WriteLine("nop");
            sw.WriteLine("pop {r4-r7}");
            sw.WriteLine("pop {pc}");
            sw.WriteLine(".endfunc");
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
            sw.WriteLine("   mov r{0},r1", len_r + 6);
            sw.WriteLine("   mov r{0},r2", len_r + 4);
            sw.WriteLine("   mov r{0},r3", len_r + 3);
            //Get Shift results
            //sw.WriteLine("   mov r3,r0");
            sw.WriteLine("  ands r1,r0");//r1=x&ShareMask2
            sw.WriteLine("  ands r3,r0");// mov r3=r0
            //Clear ALU
            sw.WriteLine("  movs r2,#0");//mov r2,#0
            //Entering the higher two shares
            for (int i = 1; i < len_r; i++)
            {
                if ((i % (shares - 1) == 0) && (i > 1))
                {
                    sw.WriteLine("   lsrs r2,r1,#4");//r2=r1>>4
                    sw.WriteLine("   lsls r1,#12");//r1=r1<<12
                    sw.WriteLine("   eors r1,r2");//r1=(r1>>4)^(r1<<12)
                    sw.WriteLine("   mov r2,r{0}", len_r + 3);//r2=AndMask
                    sw.WriteLine("   ands r1,r2");//r1=AndMask&r1
                    sw.WriteLine("   mov r{0},r1", i + 3);//r{i+3}=r1;
                }
                else
                {
                    //Right shift shares by 1=Rotated shift shares by 1
                    if (i + 3 > 7)
                    {
                        sw.WriteLine("   lsrs r2,r1,#1", i + 2);//r2=r1>>1
                        sw.WriteLine("   mov r{0},r2", i + 3);//r{i+3}=r2;
                    }
                    else
                    {
                        sw.WriteLine("   lsrs r{0},r1,#1", i + 3);//r{i+3}=r1>>1
                    }
                }
            }
            //
            sw.WriteLine("   mov r1,r{0}", len_r + 4);//ShareMask1
            sw.WriteLine("   lsrs r2,r1,#1");//r2=ShareMask1>>1
            sw.WriteLine("   eors r1,r2");//r1=ShareMask^(ShareMask1>>1)
            sw.WriteLine("  ands r1,r0");//r1=x&(ShareMask^(ShareMask1>>1))
            //Entering the lower one share
            for (int i = 1; i < len_r; i++)
            {
                if ((i % (shares - 1) == 0) && (i > 1))
                {
                    sw.WriteLine("   lsrs r2,r1,#4");//r2=r1>>4
                    sw.WriteLine("   lsls r1,#12");//r1=r1<<12
                    sw.WriteLine("   eors r1,r2");//r1=(r1>>4)^(r1<<12)
                    sw.WriteLine("   mov r2,r{0}", len_r + 3);//r2=AndMask
                    sw.WriteLine("   ands r1,r2");//r1=AndMask&r1
                    if (i + 3 > 7)
                    {
                        sw.WriteLine("   movs r0,#0");//r0=0;
                        sw.WriteLine("   mov r0,r{0}", i + 3);//r0=r{i + 3};
                        sw.WriteLine("   eors r0,r1");//r0=r0^r1;
                        sw.WriteLine("   mov r{0},r0", i + 3);//r0=r{i + 3};
                    }
                    else
                    {
                        sw.WriteLine("   eors r{0},r1", i + 3);//r0=r0^r1;
                    }
                }
                else
                {
                    sw.WriteLine("   mov r0,r{0}", len_r + 4);//ShareMask1
                    sw.WriteLine("   ands r0,r1");//r0=r1&ShareMask1
                    //Right shift shares by 1=Rotated shift shares by 1
                    if (i + 3 > 7)
                    {
                        sw.WriteLine("   lsrs r2,r0,#1");//r2=r0>>1
                        sw.WriteLine("   lsls r0,#2");//r0=r0<<2
                        sw.WriteLine("   eors r2,r0");//r2=(r1<<2)^(r1>>1)
                        sw.WriteLine("   movs r0,#0");//r0=0;
                        sw.WriteLine("   mov r0,r{0}", i + 3);//r0=r{i + 3};
                        sw.WriteLine("   eors r0,r2");//r0=r0^r2;
                        sw.WriteLine("   mov r{0},r0", i + 3);//r0=r{i + 3};
                    }
                    else
                    {
                        sw.WriteLine("   lsrs r2,r0,#1", i + 2);//r2=r0>>1
                        sw.WriteLine("   lsls r0,#2");//r0=r0<<2
                        sw.WriteLine("   eors r2,r0");//r2=(r1<<2)^(r1>>1)
                        sw.WriteLine("   eors r{0},r2", i + 3);//r{i+3}=(r1>>1)^r{i+3}
                    }
                }
            }
            //Getting the lower 1 bits

            //Sort anf as degree and delta_index
            Array.Sort(anf);
            //Use r{0} as result
            //Use r{1} as temp1
            //Use r{2} as temp2
            sw.WriteLine("   movs r0,#0");
            for (int i = 0; i < anf.Length; i++)
            {

                if (anf[i].d == 1)//linear term
                {
                    /**No Add**/
                    //if (anf[i].base_index + 3 > 7)
                    //    sw.WriteLine("   mov r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    //else
                    //{
                    //    sw.WriteLine("   mov r1,r{0}", len_r + 3);//temp1=Andmask;
                    //    sw.WriteLine("   ands r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    //}
                    sw.WriteLine("   mov r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    sw.WriteLine("   eors r0,r1");
                }
                else//only deal with degree 2 funcitons
                {
                    /**No Add**/
                    //if (anf[i].base_index + 3 > 7)
                    //    sw.WriteLine("   mov r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    //else
                    //{
                    //    sw.WriteLine("   mov r1,r{0}", len_r + 3);//temp2=Andmask;
                    //    sw.WriteLine("   ands r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    //}
                    //if (anf[i].base_index + anf[i].other_index[0] + 3 > 7)
                    //    sw.WriteLine("   mov r2,r{0}", anf[i].base_index + anf[i].other_index[0] + 3);//temp2=reg[];
                    //else
                    //{
                    //    sw.WriteLine("   mov r2,r{0}", len_r + 3);//temp2=Andmask;
                    //    sw.WriteLine("   ands r2,r{0}", anf[i].base_index + anf[i].other_index[0] + 3);//temp2=reg[];
                    //}
                    sw.WriteLine("   mov r1,r{0}", anf[i].base_index + 3);//temp1=reg[];
                    sw.WriteLine("   mov r2,r{0}", anf[i].base_index + anf[i].other_index[0] + 3);//temp2=reg[];
                    sw.WriteLine("   ands r1,r2");//temp1=temp1&temp2
                    //Xor temp to result
                    sw.WriteLine("   eors r0,r1");

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

            int slice_block = 32/size ;
            sw.WriteLine(".func TI_S_Shares");
            sw.WriteLine("TI_S_:");
            sw.WriteLine("push {lr}");
            sw.WriteLine("push {r4-r7}");
            for (int i = 0; i < 8; i++)
                sw.WriteLine("nop");
            //Pass r0 as the first share
            //Pass r1 as the second share
            //Note that and/eor/lsl/lsr do not work with high regs

            //Use r3 as result
            //Use r4 as temp1
            //Use r5 as temp2
            //Use r6 as temp3
            //Sort anf as degree and delta_index
            Array.Sort(anf);
            int sub = 0;//indicate the cycle save when using ARM instructions
            sw.WriteLine("   movs r3,#0");
            for (int i = 0; i < anf.Length; i++)
            {
                if (anf[i].d == 1)//linear term
                {
                    if (anf[i].base_index >= size)//share 1
                    {
                        sw.WriteLine("   mov r4,r1");
                        if (anf[i].base_index != size)
                        {
                            sw.WriteLine("   movs r6,#{0}",(anf[i].base_index-4)*slice_block);
                            sw.WriteLine("   rors r4,r6");
                            sub++;
                        }
                    }
                    else//share 0
                    {
                         sw.WriteLine("   mov r4,r0");
                        if(anf[i].base_index!=0)
                        {
                            sw.WriteLine("   movs r6,#{0}",(anf[i].base_index)*slice_block);
                            sw.WriteLine("   rors r4,r6");
                            sub++;
                        }
                    }
                    sw.WriteLine("   eors r3,r4");
                }
                else//only deal with degree 2 funcitons
                {
                    if (anf[i].base_index >= size)//share 1
                    {
                        sw.WriteLine("   mov r4,r1");
                        if (anf[i].base_index != size)
                        {
                            sw.WriteLine("   movs r6,#{0}", (anf[i].base_index - size) * slice_block);
                            sw.WriteLine("   rors r4,r6");
                            sub++;
                        }
                    }
                    else//share 0
                    {
                        sw.WriteLine("   mov r4,r0");
                        if (anf[i].base_index != 0)
                        {
                            sw.WriteLine("   movs r6,#{0}", (anf[i].base_index) * slice_block);
                            sw.WriteLine("   rors r4,r6");
                            sub++;
                        }
                    }
                    if ((anf[i].base_index + anf[i].other_index[0]) >= size)//share 1
                    {
                        sw.WriteLine("   mov r5,r1");
                        if ((anf[i].base_index + anf[i].other_index[0]) != size)
                        {
                            sw.WriteLine("   movs r6,#{0}", ((anf[i].base_index + anf[i].other_index[0]) - size) * slice_block);
                            sw.WriteLine("   rors r5,r6");
                            sub++;
                        }
                    }
                    else//share 0
                    {
                        sw.WriteLine("   mov r5,r0");
                        if ((anf[i].base_index + anf[i].other_index[0]) != 0)
                        {
                            sw.WriteLine("   movs r6,#{0}", (anf[i].base_index + anf[i].other_index[0]) * slice_block);
                            sw.WriteLine("   rors r5,r6");
                            sub++;
                        }
                    }
                   
                    sw.WriteLine("   ands r4,r5");//temp1=temp1&temp2
                    //Xor temp to result
                    sw.WriteLine("   eors r3,r4");

                }
            }
            sw.WriteLine("mov r0,r3");
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
            int len_e = size * (shares + 1);
            //r0 as the sharedinput
            //r1 as the highest bit version of Pmask 0'b1000 1000 1000 1000 1000
            //r2 as the highest bit mask 0'b1000 0000 0000 0000
            //r3 as temp
            //r4 as temp1
            //r5 as temp2
            sw.WriteLine("   movs r5,#0 ");
            for (int i = 0; i < shares+1; i++)
            {
                //Do Pmatrix
                sw.WriteLine("   mov r4,r1 ");
                sw.WriteLine("   mov r3,r2 ");//r3=highest bit mask
                sw.WriteLine("   ands r4,r0 ");//r4=x&r1
                //Using a full share in the previous one, so the next instruction cannot be ADD or mov rl
                sw.WriteLine("   ands r3,r4 ");//r3=x&highest bit mask;
                sw.WriteLine("   lsrs r3,#{0} ", len_e - 1 - i);//r3=r3>>(len_e-i);
                sw.WriteLine("   muls r3,r1 ");
                sw.WriteLine("   lsls r4,#4 ");//r4=r4<<4;
                sw.WriteLine("   eors r4,r3 ");//r4=r3^r4;
                sw.WriteLine("   eors r5,r4 ");//r0=r0^r4;
                if (i != shares)
                {
                    //shift r1 and r2
                    sw.WriteLine("   lsrs r1,#1 ");
                    sw.WriteLine("   lsrs r2,#1 ");
                }
            }
            sw.WriteLine("   mov r0,r5 ");
           
            sw.WriteLine("pop {pc}");
            sw.WriteLine(".endfunc");
        }
        //Write the P function in Sbox
        public void WriteASM_P_Shares(StreamWriter sw, int shares,int Pno)
        {
            int sliceblock = 32 / size;
            int P = (Pno << 1) ^ 0x01;
            int Pmask = 0;
            int t = P;
            for (int i = 0; i < size; i++)
            {
                if ((t & 0x01) == 1)
                {
                    Pmask = Pmask ^ (0x1 << (i * sliceblock));
                }
                t = t >> 1;
            }
            sw.WriteLine("Pmask={0:x}", Pmask);
            sw.WriteLine(".func TI_P_shares");
            sw.WriteLine("TI_P_shares:");
            sw.WriteLine("push {lr}");
            sw.WriteLine("push {r4-r7}");

            //r0 as the sharedinput
            //r1 as the lowest bit version of Pmask 0'x01010101
            //r2 as the highest block bits of r0 
            //r3 as the mask created by Pmask
            //r4 as temp
            //r5 as temp1
            sw.WriteLine("   movs r3,#0 ");
            sw.WriteLine("   movs r2,#1 ");
            sw.WriteLine("   lsls r2,#{0} ",sliceblock);
            sw.WriteLine("   subs r2,#1 ");
            sw.WriteLine("   lsls r2,#{0} ", 32-sliceblock);
            sw.WriteLine("   ands r2,r0 ");//get high bit
            sw.WriteLine("   lsrs r2,#{0} ", 32 - sliceblock);//get back to low bit
            sw.WriteLine("   movs r4,#1 ");
            for (int i = 0; i < sliceblock; i++)
            {
                sw.WriteLine("   movs r5,r2 ");//get high bits
                sw.WriteLine("   ands r5,r4 ");//get one high bit
                sw.WriteLine("   muls r5,r1 ");//multiply pmask
                sw.WriteLine("   eors r3,r5 ");//eor to temp mask
                sw.WriteLine("   lsls r1,#1 ");//shift pmask
                sw.WriteLine("   lsrs r2,#1 ");//shift high bits
            }
            sw.WriteLine("   lsls r0,#{0} ",sliceblock);
            
            sw.WriteLine("   eors r0,r3 ");
            sw.WriteLine("pop {r4-r7}");
            sw.WriteLine("pop {pc}");
            sw.WriteLine(".endfunc");
        }
        //Write the P function in Sbox
        public void WriteASM_P_8(StreamWriter sw, int shares,int Pno)
        {
            
            int len_e = size * (shares + 1);
            int P = (Pno<<1)^0x01;
            int Pmask = 0;
            int t = P;
            for (int i = 0; i < size; i++)
            {
                if ((t & 0x01) == 1)
                {
                    Pmask = Pmask ^ (0x8 << (i * 4));
                }
                t = t >> 1;
            }
            if(4*size<32)
                Pmask = Pmask ^ (0x8 << (size * 4));
            sw.WriteLine("Pmask={0:x}",Pmask);
            sw.WriteLine(".func TI_P");
            sw.WriteLine("TI_P:");
            sw.WriteLine("push {lr}");
            //r0 as the sharedinput
            //r1 as the highest bit version of Pmask 0'b1000 1000 1000 1000 1000
            //r2 as the highest bit mask 0'b1000 0000 0000 0000
            //r3 as the high bit mask    0'b1000 1000 1000 1000
            //r4 as temp1
            //r5 as temp2
            //r6 as temp
            sw.WriteLine("   movs r5,#0 ");
            for (int i = 0; i < shares+1; i++)
            {
                //Do Pmatrix
                sw.WriteLine("   mov r4,r3 ");
                sw.WriteLine("   mov r6,r2 ");//r6=highest bit mask
                sw.WriteLine("   ands r4,r0 ");//r4=high bits of x
                //Using a full share in the previous one, so the next instruction cannot be ADD or mov rl
                sw.WriteLine("   ands r6,r4 ");//r6=x&highest bit mask;
                sw.WriteLine("   lsrs r6,#{0} ", len_e - 1 - i);//r6=r6>>(len_e-i);
                sw.WriteLine("   muls r6,r1 ");
                sw.WriteLine("   lsls r4,#4 ");//r4=r4<<4;
                sw.WriteLine("   eors r4,r6 ");//r4=r6^r4;
                sw.WriteLine("   eors r5,r4 ");//r0=r0^r4;
                if (i != shares)
                {
                    //shift r1 and r2
                    sw.WriteLine("   lsrs r1,#1 ");
                    sw.WriteLine("   lsrs r2,#1 ");
                    sw.WriteLine("   lsrs r3,#1 ");
                }
            }
            sw.WriteLine("   mov r0,r5 ");

            sw.WriteLine("pop {pc}");
            sw.WriteLine(".endfunc");
        }
        //用ASM写出软件实现TI的代码
        public void Print_TI1b_ASM(string path, long num, string sname, int shares,byte[] TT,int Pno)
        {
            
            //Write the ASM
            String filename = String.Format(path + sname + "{0}_R{1}_{2}.S", size, round, num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            WriteASM_S_8(sw, shares, TT);
            WriteASM_P_8(sw, shares,Pno);
            sw.Close();
            fs.Close();
        }
        //用ASM写出软件实现TI的代码
        //用sharesgroup分组
        //按照相同bit进行slice，用不同寄存器存储shares
        public void Print_TI1b_ASM_SharesGroup(string path, long num, string sname, int shares, byte[] TT, int Pno)
        {

            //Write the ASM
            String filename = String.Format(path + sname + "{0}_R{1}_{2}.S", size, round, num);
            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            WriteASM_S_SharesGroup(sw, shares, TT);
            WriteASM_P_Shares(sw, shares,Pno);
            sw.Close();
            fs.Close();
        }
        public void WriteVerilog_P(StreamWriter sw, int shares, int Pno)
        {
            int P = (Pno<<1)^0x01;
            uint[,] matr=new uint[size,size];
            for(int i=0;i<size;i++)
            {
                for(int j=0;j<size;j++)
                {
                    if(j==0 && i!=size-1)
                    {
                        matr[i,j]=(uint)((((P>>(size-1-j))&0x01)>0)?1:0);
                    }
                    else
                    {
                        if(j==((i+1)%size))
                             matr[i,j]=(uint)1;
                        else
                            matr[i,j]=(uint)0;
                    }
                }
            }
            BinaryMatrix bm=new BinaryMatrix((uint)size,(uint)size);
            bm.SetValue(matr);
            bm.LeftMultiplyMatrix(bm);
            bm.GetValue(matr);
            sw.WriteLine("input[3:0] in;");
            sw.WriteLine("output[3:0] out;");
            sw.WriteLine("wire[3:0] out;");
            for(int i=0;i<size;i++)
            {
                 sw.Write("assign out[{0}]=",size-1-i);
                bool start=true;
                for(int j=0;j<size;j++)
                {
                    if(matr[i,j]==1)
                    {
                        if(start)
                            sw.Write("in[{0}]",size-1-j);
                        else
                            sw.Write("^in[{0}]",size-1-j);
                        start=false;
                    }
                    
                }
                sw.WriteLine(";");
            }
            sw.WriteLine("endmodule");
        }
        //Write the S function in Sbox
        public void WriteVerilog_S(StreamWriter sw, int shares, byte[] TT)
        {
            int[] ANF = MoebiusTrans(size * (shares - 1), TT);
            //Compute ANF terms
            ANFterm[] anf = ConvertANF(ANF, shares);
            sw.WriteLine("input[{0}:0] in;", (shares - 1) * size-1);
            sw.WriteLine("output out;");
            sw.WriteLine("wire out;");
            for (int i = 0; i < anf.Length; i++)
                sw.WriteLine("wire term_{0};", i);
            
            Array.Sort(anf);
            for (int i = 0; i < anf.Length; i++)
            {
                if (anf[i].d == 1)//linear term
                {
                    sw.WriteLine("assign term_{0}=in[{1}];", i,anf[i].base_index);
                }
                else//only deal with degree 2 funcitons
                {
                    sw.WriteLine("assign term_{0}=in[{1}]&in[{2}];", i, anf[i].base_index, anf[i].base_index+anf[i].other_index[0]);
                }
            }
            sw.Write("assign out=");
            for (int i = 0; i < anf.Length; i++)
                if (i != anf.Length - 1)
                    sw.Write("term_{0}^", i);
                else
                    sw.WriteLine("term_{0};\n", i);
            sw.WriteLine("endmodule");
        }
        //用Verilog写出硬件实现TI的代码
        //用sharesgroup分组
        public void Print_TI1b_Verilog(string path, long num, string sname, int shares, byte[] TT, int Pno)
        {

            //Write the ASM
            String filenameS = String.Format(path + sname + "S{0}_R{1}_{2}.v", size, round, num);
            FileStream fsS = new FileStream(filenameS, FileMode.Create);
            String filenameP = String.Format(path + sname + "P{0}_R{1}_{2}.v", size, round, num);
            FileStream fsP = new FileStream(filenameP, FileMode.Create);
            String Smodulename = String.Format(sname + "S{0}_R{1}_{2}", size, round, num);
            String Pmodulename = String.Format(sname + "P{0}_R{1}_{2}", size, round, num);
            StreamWriter swS = new StreamWriter(fsS);
            StreamWriter swP = new StreamWriter(fsP);
            swS.WriteLine("module  {0}(in,out);", Smodulename);
            WriteVerilog_S(swS, shares, TT);
            swP.WriteLine("module  {0}(in,out);", Pmodulename);
            WriteVerilog_P(swP, shares, Pno);
            //WriteASM_P_Shares(sw, shares);
            swS.Close();
            fsS.Close();
            swP.Close();
            fsP.Close();
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
            swScript.WriteLine("compile_ultra");
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
        //写出高位地址被AES Sbox obfuscated 的Large shares Sbox
        public void WriteObfuscatedTable_SharedGroup(int[] TITT, string outfile, int shares)
        {
            int full = shares * size;
            int[] sbox = {
  //0     1    2      3     4    5     6     7      8    9     A      B    C     D     E     F
  0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
  0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
  0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
  0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
  0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
  0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
  0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
  0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
  0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
  0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
  0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
  0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
  0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
  0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
  0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
  0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 };
            FileStream fs = new FileStream(outfile, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            
            int[] SharedSbox = new int[(0x1 << full)];
            sw.WriteLine("uint_8 SharedSbox[{0}]={{", SharedSbox.Length);
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
                        for (int m = 0; m < shares - 1; m++)
                        {
                            t = t ^ sharedinput[(m + j + shares - 2) % shares];
                            if (m != shares - 2)
                                t = t << size;
                        }
                        //Compute output
                        if (TITT[t] == 1)
                            sharedoutput[j] = SetBit(sharedoutput[j], k);
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
                //input=i,output=tempv
                int high_addr = (0xff0&i)>>4;
                int new_high_addr = sbox[high_addr];
                SharedSbox[(new_high_addr << 4) ^ (i & 0xf)] = tempv;
            }
            for (int i = 0; i < (0x1 << full); i++)
            {
                if (i == (0x1 << full) - 1)
                    sw.WriteLine("{0}}};", SharedSbox[i]);
                else
                {
                    if (i % 16 == 15)
                        sw.WriteLine("{0},", SharedSbox[i]);
                    else
                        sw.Write("{0},\t", SharedSbox[i]);
                }
            }
            sw.Close();
            fs.Close();
            
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
        //SITI+*4P
        public void OneRoundTrans_SITIM4P_Last(int[] table, int[] S)
        {

            for (int i = 0; i < len; i++)
            {
                table[i] = S[table[i]];
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
            Parallel.For(0, length, new ParallelOptions { MaxDegreeOfParallelism = 1 }, num =>
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

                for (int i = 0; i < round-1; i++)
                {
                    OneRoundTrans_SITIM4P(table, Stable[Sno], Pno);
                }
                //Remove Last P
                //OneRoundTrans_SITIM4P_Last(table, Stable[Sno]);
                //Keep Last P
                OneRoundTrans_SITIM4P(table, Stable[Sno], Pno);
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
                        sw.WriteLine("Full S Table:");
                        for (int j = 0; j < (0x1 << size); j++)
                            sw.Write("{0:x2}\t", table[j]);
                        sw.WriteLine("\n");
                        sw.WriteLine("Num of S={0}", num);

                        //Get TI bit
                        int[] ANF = MoebiusTrans(size, TT[Sno]);//Get ANF
                        //Compute the TI form from ANF
                        TI_wrapper tw = new TI_wrapper(size, shares, ANF);
                        tw.Compute_From_ANF();
                       
                        int[] TITT = tw.Get_TruthTable_F1_SharesGroup();
                        if (!VerifyTITable_SharesGroup(TITT, Stable[Sno], shares))
                        {
                            System.Console.WriteLine("TI Sbox results Incorrect!");
                            return;
                        }

                        byte[] TITTb=TTTransformation(TITT);
                        
                         //Write ASM TI implementation
                        Print_TI1b_ASM_SharesGroup(path, num, sname, 3, TITTb,Pno);
                        Print_TI1b_Verilog(path, num, sname, 3, TITTb, Pno);
                       
                        //评估代价
                        WriteSIM4PSoftwareCost_TI1b(sw, TITTb, Pno);
                        WriteSIM4PCost_Unprotected(sw, TT[Sno], Pno);
                        sw.WriteLine("******************");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.WriteLine("\n");
                        sw.Flush();
                        //写出Sbox
                        Print_FullSbox_Unprotected(path, table, num,sname);
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
