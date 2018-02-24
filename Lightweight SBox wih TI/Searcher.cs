using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Lightweight_SBox_wih_TI
{
    class Searcher
    {
        public int size;//输入长度
        public int round;//迭代轮数
        public int len;//表长
        FileStream fsScript;
        StreamWriter swScript;
        public Searcher(int s,int r)
        {
            size = s;
            round = r;
            len = (int)Math.Pow(2, size);
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
        /// 取比特函数：取input的第pos位，其中第0位处于最右边，即LSB
        /// </summary>
        /// <param name="input"></param>
        /// <param name="pos"></param>
        /// <returns></returns>
        public int GetBit(int input, int pos)
        {
            return ((input >> pos) & 0x1);
        }
        /// <summary>
        /// 置比特函数：取input的第pos位，其中第0位处于最右边，即LSB
        /// </summary>
        /// <param name="input"></param>
        /// <param name="pos"></param>
        /// <returns></returns>
        public int SetBit(int input, int pos)
        {
            return (input | (0x1<<pos));
        }
        /// <summary>
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

        /// <summary>
        /// 计算S盒的差分均匀度，sbox是其真值表表示，size是其比特级规模
        /// </summary>
        /// <param name="sbox"></param>
        /// <param name="size"></param>
        /// <returns></returns>
        public int DiffUniform(int[] sbox, int size)
        {
            int maxDiff;

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
        public int Nonlinear(int[] sbox)
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
        /// 一轮变换对应的真值表更新：
        /// size为S盒的大小
        /// table为真值表，大小为2^size
        /// func为布尔函数的真值表
        /// </summary>
        /// <param name="size"></param>
        /// <param name="table"></param>
        /// <param name="func"></param>
        public void OneRoundTrans(int[] table, int func)
        {
            int newX0, newX1, newX2, newX3;

            for (int i = 0; i <len; i++)
            {
                newX3 = GetBit(func, (table[i]) >> 1) ^ GetBit(table[i], 0);
                newX2 = GetBit(table[i], 3);
                newX1 = GetBit(table[i], 2);
                newX0 = GetBit(table[i], 1);

                table[i] = ((newX3 << 3) | (newX2 << 2) | (newX1 << 1) | newX0);
            }
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
        public int InvMoebiusTrans(int varNum, int[] ANF)
        {
            int[] ANF1=new int[ANF.Length];
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
            int func = 0;
            for (int i = 0; i < Power(2, varNum); i++)
            {
                if(ANF1[i]==1)
                    func = SetBit(func, i);
            }
            return func;
        }
        /// <summary>
        /// 寻找最优S盒：
        /// table表示S盒的真值表
        /// </summary>
        /// <param name="table"></param>
        /// <param name="size"></param>
        /// <param name="round"></param>
        /// <param name="optimalDiff"></param>
        /// <param name="optimalNonlinear"></param>
        /// <returns></returns>
        public int[] SearchOptimal(string filename, string scriptfilename,int optimalDiff, int optimalNonlinear)
        {
            int[] table = new int[len];
            int[] errorTable = new int[len];
            int count = 0;
            int[] item = new int[round];

            fsScript = new FileStream(scriptfilename, FileMode.Create);
            swScript = new StreamWriter(fsScript);

            //写Script头
            swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            swScript.WriteLine("set synthetic_library generic.sdb");                   
            swScript.WriteLine("set target_library \"fast.db slow.db\""); 
            swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");
            
            

            FileStream fs = new FileStream(filename, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);


            //每轮的布尔函数都是一样的
            int varNum = 3;
            int[] ANF = new int[Power(2, varNum)];

            for (int num = 0; num < 256; num++)
            {
                MoebiusTrans(varNum, num, ANF);

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
                    OneRoundTrans(table, item[i]);
                }

                //求解最终S盒的差分均匀度、非线性度
                int diff = DiffUniform(table, 4);
                int nonLinear = Nonlinear(table);

                if (diff == optimalDiff && nonLinear == optimalNonlinear)
                {
                    System.Console.WriteLine("diff={0},nL={1},num={2}", diff, nonLinear, num);
                    count++;
                    sw.WriteLine("{0}",num);

                    swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/F_I3O1_D2_" + num + ".v\"}");
                    swScript.WriteLine("compile -exact_map");
                    swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/areareport_" + num + ".txt\"");
                    swScript.WriteLine("remove_design -designs");

                    //return table;
                }
            }
            System.Console.WriteLine("共有{0}个最优解",count);
            sw.Close();
            fs.Close();

            swScript.WriteLine("quit");
            swScript.Close();
            fsScript.Close();

            return errorTable;
            //释放item
        }
        //从DC的areareport中提取areacost
        public double ReadAreaCost(string filename)
        {
            FileStream fs = new FileStream(filename, FileMode.Open);
            StreamReader sr = new StreamReader(fs);
            string line = sr.ReadLine();
            double cost=0;
            while (line != null)
            {
                if (line.StartsWith("Total cell area:"))
                {
                    //int e=line.LastIndexOf(' ');
                    int s=line.Substring(0,line.Length-1).LastIndexOf(' ')+1;
                    cost = Double.Parse(line.Substring(s));
                };
                line = sr.ReadLine();
            }
            sr.Close();
            fs.Close();
            return cost;
        }
        //最后一个"_"和"."之间的部分为num
        public void GetCost_GE(string path,string CostFile)
        {
            FileStream fs = new FileStream(CostFile, FileMode.Create);
            StreamWriter sw = new StreamWriter(fs);
            DirectoryInfo TheFolder = new DirectoryInfo(path);
            //遍历文件
            foreach (FileInfo NextFile in TheFolder.GetFiles())
            {
                double cost=ReadAreaCost(NextFile.FullName)/9.9792;
                string fname = NextFile.Name;
                int s = fname.LastIndexOf('_')+1;
                int e = fname.LastIndexOf('.');
                int num=Int32.Parse(fname.Substring(s,e-s));
                sw.WriteLine("{0}\t{1}", num, cost);
            }
            sw.Close();
            fs.Close();
        }
    }
}
