import java.io.File;
import java.util.Arrays;
import java.util.Scanner;
@Grab(group='org.apache.commons', module='commons-math3', version='3.6.1')
import org.apache.commons.math3.stat.inference.ChiSquareTest;

// input.vcf(.gz) mother father
public class FilterVCFforMap {

	public static void main(String[] args) throws Exception {
		Scanner s = new Scanner(new File(args[0]));
		String m = args[1]; //mother
		String f = args[2]; //father
		double p = Double.parseDouble(args[3]);
		ChiSquareTest chi = new ChiSquareTest();
		int m_index = 0, f_index = 0;
		System.err.print("CHROM\t"+"POS\t"+"calls\t"+"no_calls\t"+"AA\t"+"Aa\t"+"aa\t"+"chi_val\t"+"p_val\t"+"skip\n");
		while(s.hasNext()) {
			boolean skip=false;
			String line = s.nextLine();
			//Comment line, just print it
			if(line.startsWith("##")) System.out.print(line+"\n");
			//Header line, get parents indexes and print it
			else if(line.startsWith("#")) {
				String[] data = line.split("\t");
				m_index = Arrays.asList(data).indexOf(m);
				f_index = Arrays.asList(data).indexOf(f);
				System.out.print(line+"\n");
			}
			//SNP line, process it
			else {
				String[] data = line.split("\t");
				String CHROM = data[0];
				String POS = data[1];
				String REF = data[3];
				String ALT = data[4];
				int gt_index = Arrays.asList(data[8].split(":")).indexOf("GT");
				//Get parents genotypes
				String gt_m = data[m_index].split(":")[gt_index];
				String gt_f = data[f_index].split(":")[gt_index];
				//If one or two of the parents missing, skip the SNP
				if(gt_m.equalsIgnoreCase("./.") || gt_f.equalsIgnoreCase("./.")) skip=true;
				//If both of the parents are hom, skip the SNP\
				if((gt_m.equalsIgnoreCase("0/0")&&gt_f.equalsIgnoreCase("0/0")) ||
				   (gt_m.equalsIgnoreCase("1/1")&&gt_f.equalsIgnoreCase("1/1")) ||
				   (gt_m.equalsIgnoreCase("1/1")&&gt_f.equalsIgnoreCase("0/0")) ||
				   (gt_m.equalsIgnoreCase("0/0")&&gt_f.equalsIgnoreCase("1/1"))
				   ) skip=true;
				//Check progeny genotypes inconsistencies
				for(int i=9; i<data.length; i++) {
					if(i==m_index || i==f_index) continue;
					String gt = data[i].split(":")[gt_index];
					if(gt.equalsIgnoreCase("0/0") && (
							(gt_m.equalsIgnoreCase("1/1")&&gt_f.equalsIgnoreCase("0/1")) ||
							(gt_m.equalsIgnoreCase("0/1")&&gt_f.equalsIgnoreCase("1/1")) )
							) {
						data[i] = "./.";
					}
					else if(gt.equalsIgnoreCase("1/1") && (
							(gt_m.equalsIgnoreCase("0/0")&&gt_f.equalsIgnoreCase("0/1")) ||
							(gt_m.equalsIgnoreCase("0/1")&&gt_f.equalsIgnoreCase("0/0")) )
							) {
						data[i] = "./.";
					}
				}
				//Mendelian segregation distortion
				int calls=0, no_calls=0, AA=0, Aa=0, aa=0;
				double chi_val=0.0, p_val=0.0;
				for(int i=9; i<data.length; i++) {
					if(i==m_index || i==f_index) continue;
					String gt = data[i].split(":")[gt_index];
					if(gt.equalsIgnoreCase("0/1")) Aa++;
					else if(gt.equalsIgnoreCase("0/0")&&gt_m.equalsIgnoreCase("0/0")&&gt_f.equalsIgnoreCase("0/1")) AA++;
					else if(gt.equalsIgnoreCase("1/1")&&gt_m.equalsIgnoreCase("1/1")&&gt_f.equalsIgnoreCase("0/1")) AA++;
					else if(gt.equalsIgnoreCase("0/0")&&gt_m.equalsIgnoreCase("0/1")&&gt_f.equalsIgnoreCase("0/0")) aa++;
					else if(gt.equalsIgnoreCase("1/1")&&gt_m.equalsIgnoreCase("0/1")&&gt_f.equalsIgnoreCase("1/1")) aa++;
					else if(gt.equalsIgnoreCase("0/0")&&gt_m.equalsIgnoreCase("0/1")&&gt_f.equalsIgnoreCase("0/1")) AA++;
					else if(gt.equalsIgnoreCase("1/1")&&gt_m.equalsIgnoreCase("0/1")&&gt_f.equalsIgnoreCase("0/1")) aa++;
					else no_calls++;
				}
				calls = AA + Aa + aa;
				//System.err.print(CHROM+"\t"+POS+"\t"+calls+"\t"+no_calls+"\t"+AA+"\t"+Aa+"\t"+aa+"\t"+chi_val+"\t"+p_val+"\n");
				if(calls>0) {
					//hkxhk, F2 AA-Aa-aa
					if(gt_m.equalsIgnoreCase("0/1")&&gt_f.equalsIgnoreCase("0/1")) {
						chi_val = chi.chiSquare(new double[] {0.25*calls,0.5*calls,0.25*calls}, new long[] {AA,Aa,aa});
						p_val = chi.chiSquareTest(new double[] {0.25*calls,0.5*calls,0.25*calls}, new long[] {AA,Aa,aa});
					}
					//nnxnp, male backcross AA-Aa-0
					else if(AA!=0 && Aa!=0 && aa==0) {
						chi_val = chi.chiSquare(new double[] {0.5*calls,0.5*calls}, new long[] {AA,Aa});
						p_val = chi.chiSquareTest(new double[] {0.5*calls,0.5*calls}, new long[] {AA,Aa});
					}
					//lmxll, female backcross 0-Aa-aa
					else if(AA==0 && Aa!=0 && aa!=0) {
						chi_val = chi.chiSquare(new double[] {0.5*calls,0.5*calls}, new long[] {Aa,aa});
						p_val = chi.chiSquareTest(new double[] {0.5*calls,0.5*calls}, new long[] {Aa,aa});
					}
					//0-Aa-0
					else if(AA==0 && Aa!=0 && aa==0) {
						p_val = 1.0;
						skip=true;
					}
					//AA-0-aa
					else if(AA!=0 && Aa==0 && aa!=0) {
						p_val = 0.0;
						skip=true;
					}
				}
				else {
					skip=true;
				}
				System.err.print(CHROM+"\t"+POS+"\t"+calls+"\t"+no_calls+"\t"+AA+"\t"+Aa+"\t"+aa+"\t"+chi_val+"\t"+p_val+"\t"+skip+"\n");
				
				//Print the filtered SNP
				if(!skip && !Double.isNaN(p_val) && Double.isFinite(p_val) && p_val>=p) {
					System.out.print(String.join("\t", data)+"\n");
				}
			}
		}
		s.close();
	}
}

