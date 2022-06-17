import numpy as np
import csv
import mpmath as mp
import os 
import sys, getopt

mp.prec = 336 # binary precision 
mp.dps = 250 # decimal precision 

def Cauchy_p_calculation(input_name, output_name):
    my_data = np.genfromtxt(input_name, delimiter=',')
    p_values = my_data[1:]
    output = []
    T_list = []
    for i in range(len(p_values)):
        pval = p_values[i]
        pval = pval[np.logical_not(np.isnan(pval))]
        num_pvalue = len(pval)
        T1 = [(mp.fdiv(mp.mpf(1),mp.mpf(num_pvalue)))*mp.tan((mp.fsub(mp.mpf('0.5'), mp.mpf(p)))*mp.pi) for p in pval]
        T = sum(T1)
        p_cauchy = mp.fsub(mp.mpf('0.5'),(mp.atan(T)/mp.pi))
        p_cauchy
        output.append(p_cauchy)
        T_list.append(T)
        rounded = [float(mp.nstr(p, 5)) for p in output]
        

        with open(output_name, "w+", newline = '') as file: 
            writer = csv.writer(file)
            for (element1, element2) in zip(T_list, rounded):
                writer.writerow([element1, element2])

# read from directory and compute Cauchy p-values 


def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print('Cauchy_pval.py -i <input_path> -o <output_path>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('Cauchy_pval.py -i <input_path> -o <output_path>')
         sys.exit()
      elif opt in ("-i", "--i"):
         input_path = arg
      elif opt in ("-o", "--o"):
         output_path = arg
   print('Input path is "', input_path)
   print('Output path is "', output_path)
   
   file_names = os.listdir(input_path)
   
   for file in file_names:
     Cauchy_p_calculation(input_name = input_path+file, output_name = output_path+file)

if __name__ == "__main__":
   main(sys.argv[1:])

