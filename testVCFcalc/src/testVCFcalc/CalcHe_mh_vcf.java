package testVCFcalc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.samtools.util.CloseableIterator;

public class CalcHe_mh_vcf {
	
	// calculate expected heterozygosity for a microhap
	// uses only genotypes for inds that are phased and non-missing
	// @param vcList contains the VariantContexts of SNPs
	//   phase with WhatsHap that are in the window
	public static double[] calcHe(ArrayList<VariantContext> vcList) {
		double[] ret = new double[]{-1, 0, 0, 0}; // expected Het, number of individuals used, num failed geno, num invalid phase
		if (vcList.size() < 1) return ret; // make sure there are genotypes
		// first, initiate genotypes
		ArrayList<HaplotypeG> genos = new ArrayList<HaplotypeG>(vcList.get(0).getNSamples());
		for(int i = 0; i < vcList.get(0).getNSamples(); i++) genos.add(new HaplotypeG(vcList.size()));
		
		// now build up the genotypes from the constituent SNPs
		for(VariantContext vc : vcList) {
			int i = -1; // index of individual in genos
			for(Genotype g : vc.getGenotypesOrderedByName()) {
				i++;
				if(!genos.get(i).isCalled) continue;
				if(g.isNoCall()) {
					genos.get(i).setMissing(); // genotype is no longer good to be included
					continue;
				}
				
				// if we're good, then add to haplotypes
				genos.get(i).a1.add(g.getAllele(0).getBases()[0]);
				genos.get(i).a2.add(g.getAllele(1).getBases()[0]);
				if(g.isHet()) genos.get(i).phaseSet.add(g.isPhased() ? g.getExtendedAttribute("PS") : ".");	
			}	
		}
		// calculate He
		Map<List<Byte>, Integer> af = new HashMap<List<Byte>, Integer>();
		for(HaplotypeG g : genos) {
			// count if haplotypes are called and fully phased
			if(!g.isCalled) {
				ret[2]++;
			} else if(g.validPhase()){
				ret[1]++; // samples with useable genotypes
				af.put(g.a1, af.getOrDefault(g.a1, 0) + 1); // tabulate allele frequencies
				af.put(g.a2, af.getOrDefault(g.a2, 0) + 1);
			} else {
				ret[3]++; // called but unphased
			}
		}
		if(ret[1] == 0) return ret; // if no useable genotypes, can't estimate
		double nAlleles = 2 * ret[1];
		ret[0] = 1;
		for(int v : af.values()) {
			ret[0] -= Math.pow(v / nAlleles, 2);
		}
		return ret;
	}
	
	

	public static void main(String[] args) throws IOException {
		int wSinput = 125;
		String vcfInPath = null;
		int flag = 0;
		while(flag < args.length) {
			if(args[flag].equals("-v")) {
				vcfInPath = args[++flag];
			} else if (args[flag].equals("-w")) {
				wSinput = Integer.parseInt(args[++flag]);
			} else {
				System.out.println("Error: command line argument: " + 
						args[flag] + " is not recognized.");
				return;
			}
			flag++;
		}
		final int wS = wSinput; // making final to pass to predicate

		System.out.println("vcf file Input: " + vcfInPath);
		System.out.println("window size: " + wS);
		
		// Open VCF
		final File vcfFile = new File(vcfInPath);
		// VCFFileReader(java.io.File file, boolean requireIndex)
		VCFFileReader vcf = new VCFFileReader(vcfFile, false);		
		CloseableIterator <VariantContext> vcfIter = vcf.iterator();
		VariantContext vc = vcfIter.next(); // get first variant
		ArrayList<VariantContext> curSNP = new ArrayList<VariantContext>(); // SNP positions in current window
		while(!vc.getType().equals(VariantContext.Type.SNP)  && vcfIter.hasNext()) {
			vc = vcfIter.next();
		}
		if(vc.getType().equals(VariantContext.Type.SNP)) curSNP.add(vc); // using if in case no SNPs in the whole file
		
//		int testingCounter = 0;
		BufferedWriter fileOut = Files.newBufferedWriter(Paths.get("vcf_he_output.txt"));
		fileOut.write("Chr\tPos\tHe\tnumInds\tnumFailGeno\tnumInvalidPhase");
		fileOut.newLine();
		while(vcfIter.hasNext()) {
			vc = vcfIter.next();
			
			// determine if next SNP is in window
			if(curSNP.get(0).getContig().equals(vc.getContig()) && (vc.getStart() - curSNP.get(0).getStart()) < wS) {
				// add to window if the variant is a SNP
				if(vc.getType().equals(VariantContext.Type.SNP)) curSNP.add(vc);
			} else {
				// estimate heterozygosity
				double[] res = calcHe(curSNP);
				// build up line of output
				StringBuilder line = new StringBuilder(300); // 300 is too much, but really not THAT much memory
				// chr pos He numInds numFail numInvalidPhase
				line.append(curSNP.get(0).getContig());
				line.append("\t");
				for (int i=0, m = curSNP.size() - 1; i < curSNP.size(); i++) {
					line.append(curSNP.get(i).getStart());
					if(i < m) line.append(",");
				}
				line.append("\t");
				line.append(res[0]);
				line.append("\t");
				line.append(res[1]);
				line.append("\t");
				line.append(res[2]);
				line.append("\t");
				line.append(res[3]);
				// now write out
				fileOut.write(line.toString());
				fileOut.newLine(); // add newline
				// start a new window with at least one new SNP
				if(curSNP.get(0).getContig().equals(vc.getContig())) {
					// same chromosome, so preserve SNPs within new window
					final int newStart = vc.getStart();
					curSNP.removeIf(x -> ((newStart - x.getStart()) >= wS));
				} else {
					// diff chromosome, start all over
					curSNP.clear();
				}
				curSNP.add(vc); // add the new SNP

			}
			
//			testingCounter++;
//			if(testingCounter > 300) break; 	
		}
		
		// handle the last window
		
		// estimate heterozygosity
		double[] res = calcHe(curSNP);
		// build up line of output
		StringBuilder line = new StringBuilder(300); // 300 is too much, but really not THAT much memory
		// chr pos He numInds numFail numInvalidPhase
		line.append(curSNP.get(0).getContig());
		line.append("\t");
		for (int i=0, m = curSNP.size() - 1; i < curSNP.size(); i++) {
			line.append(curSNP.get(i).getStart());
			if(i < m) line.append(",");
		}
		line.append("\t");
		line.append(res[0]);
		line.append("\t");
		line.append(res[1]);
		line.append("\t");
		line.append(res[2]);
		line.append("\t");
		line.append(res[3]);
		// now write out
		fileOut.write(line.toString());
		fileOut.newLine();
		
		fileOut.close();
		vcfIter.close();
		vcf.close();
	}

}
