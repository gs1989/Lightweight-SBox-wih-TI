using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Lightweight_SBox_wih_TI
{
    class TI_Searcher//搜索d个变元的乘积齐次项的TI
    {
        public int d;//搜索项的次数,即变元个数
        public int s;//share 个数(s>=d+1)

        public int[] coef;//系数，共有(s)^d个项，为1时表示有效(f0的系数,系数表示这个变元的下标)

        public TI_Searcher(int di, int si)
        {
            if (si < di + 1)
            {
                System.Console.WriteLine("参数有误,次数={0}, share个数={1}", di, si);
                si = di;
            }
            d = di;
            s = si;
            int num = (int)Math.Pow(s, d);
            coef = new int[num];
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
        //将一个d长的s进制数转换为整数
        public int Trans_Back(int[] t, int s)
        {
            int result = 0;
            for (int i = t.Length - 1; i > -1; i--)
            {
                result = result + t[i];
                if (i != 0)
                    result *= s;
            }
            return result;
        }
        public void ComputeCoef()
        {
            int delta_num = (int)Math.Pow(s, d - 1);
            int[] index = new int[d];
            index[0] = 1;
            for (int delta = 0; delta < delta_num; delta++)
            {
                int[] t = Trans(delta, s, d - 1);
                bool[] Used = new bool[s];
                Used[1] = true;
                for (int k = 1; k < d; k++)
                {
                    index[k] = (index[0] + t[k - 1]) % s;
                    Used[index[k]] = true;
                }
                //找一个没有用到的下标
                int em = 0;
                while (Used[em])
                    em++;
                //计算 (em+m)%s=0;
                if (em == 0)
                    em = s;
                //此时+(s-em)即为所求，缺少0下标的项
                for (int k = 0; k < d; k++)
                {
                    index[k] = (index[k] + s - em) % s;
                }
                //转换回去
                int ind = Trans_Back(index, s);
                //存入coef数组
                if (coef[ind] != 0)
                {
                    System.Console.WriteLine("错误！不同轮间出现重复项！");
                    return;
                }
                coef[ind] = 1;
            }

        }
        //获取正常v个变元时，对应的coef数组
        //输入v为变元数;ContainVars为v长的数组，其中为1的位置表示当前Searcher中有的变元，为0的表示没有的变元
        //返回的coef数组元素取值范围为[0,s]，0表示没有这个变元
        public int[] GetCorrectTerm(int[] ContainVars, int v)
        {

            int[] newcoef = new int[(int)Math.Pow(s + 1, v)];

            int[] Corres = new int[d];
            int ind = 0;
            for (int i = 0; i < v; i++)
                if (ContainVars[i] != 0)
                {
                    Corres[ind] = i;//表示当前类中序号为ind的变元对应原来第i个变元
                    ind++;
                }
            for (int k = 0; k < coef.Length; k++)
            {
                if (coef[k] == 1)
                {
                    int[] temp = Trans(k, s, d);
                    int[] temp_ori = new int[v];
                    for (int i = 0; i < d; i++)//增添有的项
                    {
                        temp_ori[Corres[i]] = temp[i] + 1;
                    }
                    int result = Trans_Back(temp_ori, s + 1);//转换回去，注意此时0表示没有变元!!!!!
                    newcoef[result] = 1;
                }
            }
            return newcoef;
        }

        //循环移位1
        public void Shift1()
        {
            int[] newcoef = new int[(int)Math.Pow(s, d)];

            for (int k = 0; k < coef.Length; k++)
            {
                if (coef[k] == 1)
                {
                    int[] temp = Trans(k, s, d);
                    for (int i = 0; i < d; i++)//增添有的项
                    {
                        temp[i] = (temp[i] + 1) % s;
                    }
                    int result = Trans_Back(temp, s);//转换回去，注意此时0表示没有变元!!!!!
                    newcoef[result] = 1;
                }
            }
            for (int k = 0; k < coef.Length; k++)
                coef[k] = newcoef[k];
        }




    }
}
