package testVCFcalc;

import java.util.ArrayList;
import java.util.List;
import java.util.HashSet;
import java.util.Set;

public class HaplotypeG {
	// class holding (micro)haplotype values for diploids
	// not particularly refined, just quick to write
	List<Byte> a1;
	List<Byte> a2;
	boolean isCalled; // false to indicate missing genotype missing genotype
	Set <Object> phaseSet; // PS tags
	int numHet = 0; // number of heterozygous positions (for detecting if phased or not)
	
	public HaplotypeG() {
		this.a1 = new ArrayList<Byte>();
		this.a2 = new ArrayList<Byte>();
		this.isCalled = true;
		this.phaseSet = new HashSet<Object>();
	}
	public HaplotypeG(int initialStringSize) {
		this.a1 = new ArrayList<Byte>(initialStringSize);
		this.a2 = new ArrayList<Byte>(initialStringSize);
		this.isCalled = true;
		this.phaseSet = new HashSet<Object>();
	}
	public HaplotypeG(int initialStringSize, boolean initialIsCalled) {
		this.a1 = new ArrayList<Byte>(initialStringSize);
		this.a2 = new ArrayList<Byte>(initialStringSize);
		this.isCalled = initialIsCalled;
		this.phaseSet = new HashSet<Object>();
	}
	
	public boolean isHom() {
		return this.a1.equals(this.a2);
	}
	
	public void setMissing() {
		this.isCalled = false;
		this.a1.clear();
		this.a2.clear();
		this.phaseSet.clear();
	}
	
	public boolean validPhase() {
		if(this.phaseSet.size() > 1) return false; // more than one phase set
		// unphased het genotypes
		if(this.numHet > 1 && (this.phaseSet.size() == 0 || this.phaseSet.contains("."))) return false;
		return true;
	}
	
	public int size() {
		return this.a1.size();
	}
	
}
