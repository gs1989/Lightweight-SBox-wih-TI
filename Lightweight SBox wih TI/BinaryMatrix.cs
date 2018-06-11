using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Lightweight_SBox_wih_TI
{
    class BinaryMatrix
    {
        private uint[,] value;
        private uint rows;
        private uint columns;
        /***************************/
        //以下是基本操作
        /***************************/
        //构造方法:申请空间
        public BinaryMatrix(uint r, uint c)
        {
            rows = r;
            columns = c;
            value = new uint[r, c];
        }
        //取行数
        public uint getRowsNum()
        {
            return rows;
        }
        //取列数
        public uint getColumnNum()
        {
            return columns;
        }
        //获取内容
        public void GetValue(uint[,] t)
        {
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    t[i, j] = value[i, j];
        }
        //设置内容
        public void SetValue(uint[,] t)
        {
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    value[i, j] = t[i, j];
        }
        //重置矩阵
        public void ClearMatrix()
        {
            rows = 0;
            columns = 0;
            value = null;
        }
        //重新申请空间
        public void Reinit(uint r, uint c)
        {
            rows = r;
            columns = c;
            value = new uint[r, c];
        }
        /***************************/
        //以下是矩阵运算
        /***************************/
        //左乘二元向量
        public void LeftMultiplyVec(uint[] t)
        {
            //检查长度是否符合要求
            if (t.Length != columns)
            {
                Console.WriteLine("Length error! tlen={0},columns number={1}", t.Length, columns);
            }
            uint[] temp = new uint[rows];
            for (int i = 0; i < rows; i++)
            {
                temp[i] = 0;
                for (int j = 0; j < columns; j++)
                    temp[i] ^= value[i, j] * t[j];
            }
            Array.Copy(temp, t, columns);
        }
        //左乘二元向量
        public uint[] LeftMultiplyVec_Return(uint[] t)
        {
            //检查长度是否符合要求
            if (t.Length != columns)
            {
                Console.WriteLine("Length error! tlen={0},columns number={1}", t.Length, columns);
            }
            uint[] temp = new uint[rows];
            for (int i = 0; i < rows; i++)
            {
                temp[i] = 0;
                for (int j = 0; j < columns; j++)
                    temp[i] ^= value[i, j] * t[j];
            }
            return temp;
        }
        //左乘矩阵
        public void LeftMultiplyMatrix(BinaryMatrix t1)
        {
            uint[,] t = new uint[t1.getRowsNum(), t1.getColumnNum()];
            t1.GetValue(t);
            //检查长度是否符合要求
            if (t1.getRowsNum() != columns)
            {
                Console.WriteLine("Length error! t1={0}*{1},this={2}*{3}", t1.getRowsNum(), t1.getColumnNum(), getRowsNum(), getColumnNum());
            }
            uint[,] temp = new uint[getRowsNum(), t1.getColumnNum()];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < t1.getColumnNum(); j++)
                {
                    temp[i, j] = 0;
                    for (int k = 0; k < columns; k++)
                        temp[i, j] ^= value[i, k] * t[k, j];
                }
            }
            uint r = rows;
            uint c = t1.getColumnNum();
            t1.ClearMatrix();
            t1.Reinit(r,c);
            t1.SetValue(temp);

        }
    }
}
