using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
namespace Lightweight_SBox_wih_TI
{
    class Program
    {
        static int[] GetTruthTable(long num, int n)
        {
            long len = (long)Math.Pow(2, n);
            int[] truth = new int[len];
            for (int i = 0; i < len; i++)
                truth[i] = GetBit(num, i);
            return truth;
        }
       
        static int GetBit(long x, int no)
        {
            long mask = 0x01;
            mask = mask << no;
            if ((x & mask) != 0)
                return 1;
            else
                return 0;
        }
        static void Main(string[] args)
        {
            int size = 8;
            int round =5;
            //Searcher8 se = new Searcher8(size, round, 3);
            //se.SearchOptimal_ShiftInvariant_53SBitP_ShiftC("ShiftTIInvariant_n5_d2.bin", "ShiftTIInvariant_n3_d2.bin", "ShiftInvariant_53SBitP_ShiftC_R5.txt", 8, 96);
            //se.SearchOptimal("In6Out2_R6_shift3.txt", "FI6O2.sh", 8, 96);

           // se.SearchOptimal_WithP_BitPerm("In4Out4_R6_WithP_BitPerm.txt", 8, 100);
            //se.SearchOptimal_X2_BitP_FixC("X2_BitP_FixC_R10.txt", "X2_BitP_FixC_R10.sh", 8, 90);
            ShiftInvariant si = new ShiftInvariant(8, 2);//5 bit 2次函数
            si.SearchTIPermutation("ShiftTIInvariant_n8_d2_new.bin");
            //FileStream fs = new FileStream("In6Out2_R6_shift3.txt", FileMode.Open);
            //StreamReader sr = new StreamReader(fs);

            //FileStream fsScript;
            //StreamWriter swScript;

            //fsScript = new FileStream("FI6O2.sh", FileMode.Create);
            //swScript = new StreamWriter(fsScript);

            ////写Script头
            //swScript.WriteLine("set search_path \"/eda/synopsys/B-2008.09/dpaTest /eda/synopsys/B-2008.09/dw/sim_ver /eda/synopsys/B-2008.09/dw/syn_ver /eda/synopsys/B-2008.09/libraries/syn\"");
            //swScript.WriteLine("set synthetic_library generic.sdb");
            //swScript.WriteLine("set target_library \"fast.db slow.db\"");
            //swScript.WriteLine("set link_library \"fast.db slow.db dw_foundation.sldb\"");

            //String line = sr.ReadLine();
            //int[] ANFt = new int[64];
            //TI_wrapper tw = new TI_wrapper(6, 3, ANFt);
            //int shares = 3;
            //int v = 6;
            //while (line != null)
            //{
            //    if (line.Contains("内部函数真值表="))
            //    {
            //        long num = Convert.ToInt64(line.Substring(line.IndexOf('表') + 2), 16);
            //        int[] tt = GetTruthTable(num, v);
            //        int[] ANF = tw.TruthTable2ANF(tt);
            //        tw = new TI_wrapper(v, shares, ANF);
            //        tw.Compute_From_ANF();
            //        tw.Print_TruthTable_F1("F_I6O1_D2_" + Convert.ToString(num, 16) + ".v", num);

            //        swScript.WriteLine("read_file -format verilog {\"/mnt/hgfs/share/Circular/F_I6O1_D2_" + Convert.ToString(num, 16) + ".v\"}");
            //        swScript.WriteLine("compile -exact_map");
            //        swScript.WriteLine("report_area > \"/mnt/hgfs/share/Circular/areareports/areareport_" + Convert.ToString(num, 16) + ".txt\"");
            //        swScript.WriteLine("remove_design -designs");

            //    }
            //    line = sr.ReadLine();
            //}
            //sr.Close();
            //fs.Close();

            //swScript.WriteLine("quit");
            //swScript.Close();
            //fsScript.Close();

            //e.SearchOptimal("FI3O1_4GF.txt", "FI3O1_4GF.sh", 4, 4);
            //se.GetCost_GE("D:\\Projects\\Lightweight SBox wih TI\\Lightweight SBox wih TI\\bin\\Debug\\areareports\\4bit\\", "FI3O1_GECost.txt");

            //for (int i = 0; i < Power(2, size); i++)
            //{
            //    Console.Write("{0}, ", table[i]);
            //}
            //Console.WriteLine();
            //Searcher se = new Searcher(4, 4);
            //int[] Anf=new int[8];
            //se.MoebiusTrans(3, 3, Anf);
            //se.size = 0;
        }

 

    }
}
