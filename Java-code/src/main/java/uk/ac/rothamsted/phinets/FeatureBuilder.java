package uk.ac.rothamsted.phinets;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.apache.log4j.net.SMTPAppender;

import net.sourceforge.utils.DefaultSuppliers;
import uk.ac.rothamsted.phinets.util.BlastSet;
import uk.ac.rothamsted.phinets.util.DelimitedIO;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimaps;
import com.google.common.collect.SetMultimap;

public class FeatureBuilder {
	
	public static final void main(String  [] args) throws Exception{
		simplifyNetworks();
		buildAnnotation();
		domainRarityStats();
		String fa = "D:\\projects\\path_net\\seqsim\\fa\\";
		String blastFile = "D:\\projects\\path_net\\seqsim\\decypher_file\\blastP_allFungiVSallFungi_out_050516_final_01[removed_strains].txt";
		String outDir = "D:\\projects\\path_net\\seqsim\\scores\\";
		//computeBlastScores(fa, blastFile, outDir);
	}
	
	//fa = "D:\\projects\\path_net\\seqsim\\fa\\";
	//blastFile = "D:\\projects\\path_net\\seqsim\\decypher_file\\blastP_allFungiVSallFungi_out_050516_final_01.txt"
	//outDir = "D:\\projects\\path_net\\seqsim\\scores\\"
	public static final void computeBlastScores(String fastaFilesDir, String blastFile, String outDir) throws Exception{
		//String fastaFilesDir = "D:\\projects\\path_net\\seqsim\\fa\\";
		BlastSet bs = new BlastSet();
		bs.buildSubsetIndex(fastaFilesDir);
		double min = bs.loadBlast(blastFile);
		
		System.err.println("Computing scores...");
		bs.generateScores(outDir, min);
		System.err.println("Complete.");
	}
	
	public static void simplifyNetworks() throws Exception{
		String inDir = "D:\\projects\\path_net_combined\\networks\\complete\\";
		String outDir = "D:\\projects\\path_net_combined\\networks\\simple\\";

		for(File f : new File(inDir).listFiles()){
			BufferedReader br = new BufferedReader(new FileReader(f));
			String [] name =  f.getName().split("_");
			BufferedWriter bw = new BufferedWriter(new FileWriter(outDir+name[0]+"_"+name[1]));
			String line = br.readLine();
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				bw.write(data[0]+"\t"+data[1]+"\n");
			}
			br.close();
			bw.flush();
			bw.close();
		}
	}
	
	public static void domainRarityStats() throws Exception{
		String out = "D:\\projects\\path_net_combined\\features\\rarity\\";
		String domainsDir = "D:\\projects\\path_net_combined\\domains\\";
		String speciesFile = "D:\\projects\\path_net_combined\\species_list.tab";
		String anno = "D:\\projects\\path_net_combined\\features\\anno\\";
		String netDir = "D:\\projects\\path_net_combined\\networks\\complete\\";
		
		BufferedWriter bw5 = new BufferedWriter(new FileWriter(out+"total_rarity.tab"));
		bw5.write("gene\tphe   \trarity\trarity_net\timpact\timpact_adj\n");
		
		BufferedReader br1 = new BufferedReader(new FileReader(speciesFile));
		String line1;
		Map<String, String> id2ddiFile = new HashMap<String, String>();
		while ((line1 = br1.readLine()) != null) {
			String[] data1 = line1.split("\t");
			String id = data1[0];
			//System.err.println("Processing "+id);
			id2ddiFile.put(data1[0], data1[1]);
		}
		br1.close();
		double [] sum = new double[2];
		double [] counts = new double[2];
		
		double [] sum1 = new double[2];
		double [] counts1 = new double[2];
		for (Entry<String, String> ent : id2ddiFile.entrySet()) {
			//Load network
			BufferedReader br2 = new BufferedReader(new FileReader(netDir+ent.getKey()+"_interactome[complete].tab"));
			SetMultimap<String, String> net = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
			String line;
			while ((line = br2.readLine()) != null) {
				String[] data = line.split("\t");
				net.put(data[0], data[1]);
				net.put(data[1], data[0]);
			}
			br2.close();
			
			BufferedReader br = new BufferedReader(new FileReader(domainsDir+ent.getValue()));
			SetMultimap<String, String> map = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
			
			HashMultiset<String> counter = HashMultiset.create();
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				String id = data[0];
				String dom = data[1].split("\\.")[0];
				counter.add(dom);
				map.put(id, dom);
			}
			br.close();


			BufferedWriter bw = new BufferedWriter(new FileWriter(out+ent.getKey()+"_rarity.tab"));
			Map<String, String> an = DelimitedIO.readSimpleMap(anno+ent.getKey()+"_anno.tab", "\t", 0, 2);
			bw.write("gene\trarity\trarity_net\timpact\timpact_adj\n");
			
			double avgRar = 0d;
			double occs = 0d;
			for (Entry<String, Collection<String>> ent1 : map.asMap().entrySet()) {
				if(ent1.getValue().size()==1){
					avgRar += rarity(ent1.getValue(), counter, 0d);	
					occs += 1d;
				}
				
			}
			avgRar = 0d;
			//System.err.println(avgRar);
			for (String node : net.keySet()) {
				double min = rarity(map.get(node), counter, avgRar);
				double rarities = 0d;
				Collection<String> neighs = net.get(node);
				for(String neigh : neighs){
					double rr = rarity(map.get(neigh), counter, avgRar);
					rarities += rr;
				}
				if(neighs.size() == 0){
					rarities = 0d;	
				}
				else{
					rarities = rarities/((double)neighs.size());
				}
				
				double [] functionalImpact = complexityImpactTotal(node, counter, net, map, avgRar);
				
				bw.write(node+"\t"+String.valueOf(min)+"\t"+String.valueOf(rarities)+"\t"+functionalImpact[0]+"\t"+functionalImpact[1]+"\n");
				
				if(!an.containsKey(node)){
					continue;
				}
				System.err.println(node+"\t"+String.valueOf(min)+"\t"+an.get(node));
				if(an.get(node).equals("pathogenicity_related")){
					counts[0]+=1d;
					sum[0]+=functionalImpact[0];
					
					counts1[0]+=1d;
					sum1[0]+=functionalImpact[1];
					bw5.write(node+"\tpathogenicity_related\t"+String.valueOf(min)+"\t"+String.valueOf(rarities)+"\t"+functionalImpact[0]+"\t"+functionalImpact[1]+"\n");
				}
				else if(an.get(node).equals("pathogenicity_unrelated")){
					counts[1]+=1d;
					sum[1]+=functionalImpact[0];
					
					counts1[1]+=1d;
					sum1[1]+=functionalImpact[1];
					bw5.write(node+"\tpathogenicity_unrelated\t"+String.valueOf(min)+"\t"+String.valueOf(rarities)+"\t"+functionalImpact[0]+"\t"+functionalImpact[1]+"\n");
				}
			}
			bw.flush();
			bw.close();
		}
		bw5.flush();
		bw5.close();
		System.err.println("Path: "+(sum[0]/counts[0])+" Non-path: "+(sum[1]/counts[1]));
		System.err.println("Path: "+(sum1[0]/counts1[0])+" Non-path: "+(sum1[1]/counts1[1]));
	}
	
	
	public static double rarity(Collection<String> domains, HashMultiset<String> counter, double avg){
		if(domains==null){
			return avg;
		}
		int min = Integer.MAX_VALUE;
		for(String dom : domains){
			int count = counter.count(dom);
			if(count < min ){
				min = count;
			}
		}
		return 1d/((double)min);	
	}
	
	public static double [] complexityImpactTotal(String protein, HashMultiset<String> counter, SetMultimap<String, String> net, SetMultimap<String, String> doms, double avg){
		Collection<String> neighs = net.get(protein);
		double num = ((double)neighs.size())+1d;
		double [] totals = complexityImpact(doms.get(protein), counter, 1, avg);
		for(String neigh : neighs){
			double [] values = complexityImpact(doms.get(protein), counter, net.get(neigh).size(), avg);
			totals[0] += values[0];
			totals[1] += values[1];
		}
		totals[0] = totals[0]/num;
		totals[1] = totals[1]/num;
		return totals;
	}
	
	public static double [] complexityImpact(Collection<String> domains, HashMultiset<String> counter, int degree, double avg){
		double sum = 0d;
		for(String dom : domains){
			sum = sum + 1d/((double)counter.count(dom));
		}
		if(sum == 0){
			sum = avg;
		}
		return new double[]{sum, sum/((double)degree)};	
	}
	
	public static void buildAnnotation() throws Exception{
		String netDir = "D:\\projects\\path_net_combined\\networks\\complete\\";
		String phiDir = "D:\\projects\\path_net_combined\\phibase\\";
		String domainsDir = "D:\\projects\\path_net_combined\\domains\\";
		//String blastFile = "D:\\projects\\path_net\\seqsim\\decypher_file\\blastP_allFungiVSallFungi_out_050516_final_01.txt";
		String speciesFile = "D:\\projects\\path_net_combined\\species_list.tab";
		String outDir = "D:\\projects\\path_net_combined\\features\\";
		
		SetMultimap<String, String> id2innet = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		SetMultimap<String, String> id2inphi = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		Map<String, SetMultimap<String, String>> id2net = new HashMap<String, SetMultimap<String, String>>();
		HashMultiset<String> domainCount = HashMultiset.create();
		
		
		BufferedReader br1 = new BufferedReader(new FileReader(speciesFile));
		String line1;
		BufferedWriter bwStats = new BufferedWriter(new FileWriter(outDir+"statistics.tab"));
		bwStats.write("Species\tPhi-base annotations\tAnnotated genes\tNodes\tEdges\tContext-specific\tPathogenicity-related\tPathogenicity-unrelated\n");

		Map<String, String> id2ddiFile = new HashMap<String, String>();
		while ((line1 = br1.readLine()) != null) {
			String[] data1 = line1.split("\t");
			String id = data1[0];
			System.err.println("Processing "+id);
			id2ddiFile.put(data1[0], data1[1]);
			//Load network
			BufferedReader br = new BufferedReader(new FileReader(netDir+id+"_interactome[complete].tab"));
			SetMultimap<String, String> net = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
			id2net.put(id, net);
			String line;      
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				net.put(data[0], data[1]);
				net.put(data[1], data[0]);
				id2innet.put(id, data[0]);
				id2innet.put(id, data[1]);
			}
			br.close();
			
			Set<String> phis = new HashSet<String>();
			//Load annotations
			br = new BufferedReader(new FileReader(phiDir+id+"_phi.tab"));

			//maps.get("unaltered pathogenicity").removeAll(remove);
			//pathogenicity_related pathogenicity_unrelated
			//increased sensitivity to chemical -- increased resistance to chemical -- pathogen death
			Set<String> exclude = new HashSet<String>(){{add("increased sensitivity to chemical");add("increased resistance to chemical");add("pathogen death");}};
			
			
			
			SetMultimap<String, String> pheMap = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
			
			id2net.put(id, net);
			
			int annotations = 0;
			Set<String> genes = new HashSet<String>();
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				phis.add(data[1].toLowerCase());
				
				pheMap.put(data[0], data[1].toLowerCase());
				genes.add(data[0]);
				annotations++;
			}
			br.close();
			SetMultimap<String, String> phenoLists = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
			BufferedWriter bw = new BufferedWriter(new FileWriter(outDir+id+"_anno.tab"));
			bw.write("gene\tspecies\tphe\tphibase\n");
			int both = 0;
			for (Entry<String, Collection<String>> ent : pheMap.asMap().entrySet()) {
				Set<String> temp = new HashSet<String>(ent.getValue());
				String phiText = temp.toString().replaceFirst("\\[", "").replaceFirst("\\]", "");
				temp.removeAll(exclude);
				if(temp.isEmpty()){
					continue;
				}
				boolean isu = temp.remove("unaltered pathogenicity");
				String type =  null;
				if(isu && temp.isEmpty()){
					type = "pathogenicity_unrelated";
				}
				else if(isu){
					both++;
				}
				else{
					type = "pathogenicity_related";
				}
				if(type != null){
					id2inphi.put(id, ent.getKey());
					phenoLists.put(type, ent.getKey());
					bw.write(ent.getKey()+"\t"+id+"\t"+type+"\t"+phiText+"\n");
				}
			}
			
			bw.flush();
			bw.close();
			
			
			System.err.println(phis);
			String p = String.valueOf(phenoLists.get("pathogenicity_related").size());
			String u = String.valueOf(phenoLists.get("pathogenicity_unrelated").size());
			phenoLists.get("pathogenicity_related").retainAll(net.keySet());
			phenoLists.get("pathogenicity_unrelated").retainAll(net.keySet());
			String pn = String.valueOf(phenoLists.get("pathogenicity_related").size());
			String un = String.valueOf(phenoLists.get("pathogenicity_unrelated").size());
			bwStats.write(id+"\t"+String.valueOf(annotations)+"\t"+String.valueOf(genes.size())+"\t"+net.keySet().size()+"\t"+String.valueOf((net.size()/2))+"\t"+String.valueOf(both)+"\t"+p+"("+pn+")\t"+u+"("+un+")\n");
			
		}
		br1.close();
		bwStats.flush();
		bwStats.close();
		
		
		Set<String> selected = new HashSet<String>(){{
			add("Botrytis_cinerea");
			add("Magnaporthe_oryzae");
			add("Ustilago_maydis");
			add("Fusarium_graminearum");
		}};
		
		for(String sel : selected){
			SetMultimap<String, String> gene2domains = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
			BufferedReader br = new BufferedReader(new FileReader(domainsDir+id2ddiFile.get(sel)));
			String line;
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				String gene = data[0];
				String domain = data[1].split("\\.")[0];
				gene2domains.put(gene, domain);
			}
			br.close();
			
			for (Entry<String, Collection<String>> ent : gene2domains.asMap().entrySet()) {
				if(!id2inphi.get(sel).contains(ent.getKey())){
					continue;
				}
				domainCount.addAll(ent.getValue());
				for(String s : id2net.get(sel).get(ent.getKey())){
					domainCount.addAll(gene2domains.get(s));
				}
			}
		}

		BufferedWriter bw = new BufferedWriter(new FileWriter(outDir+"domain_counts.tab"));
		int num = 0;
		List<String> domainOrder = new ArrayList<String>();
		for(com.google.common.collect.Multiset.Entry<String> ent : domainCount.entrySet()){
			//948
			if(ent.getCount() < 10){
				continue;
			}
			domainOrder.add(ent.getElement()); 
			bw.write(ent.getElement()+"\t"+ent.getCount()+"\n");
			num++;
		}
		bw.flush();
		bw.close();
		
		
		
		for(String id : id2ddiFile.keySet()){
			SetMultimap<String, String> gene2domains1 = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
			BufferedReader br = new BufferedReader(new FileReader(domainsDir+id2ddiFile.get(id)));
			String line;
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				String gene = data[0];
				String domain = data[1].split("\\.")[0];
				gene2domains1.put(gene, domain);
			}
			br.close();
			SetMultimap<String, String> gene2domains2 = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
			for (Entry<String, Collection<String>> ent : gene2domains1.asMap().entrySet()) {
				String gene = ent.getKey();
				for(String s : id2net.get(id).get(gene)){
					gene2domains2.putAll(gene, gene2domains1.get(s));
				}
			}
			
			BufferedWriter bwDom1 = new BufferedWriter(new FileWriter(outDir+"domains\\"+id+"_domains_all.tab"));
			bwDom1.write("gene");
			for(String dom : domainOrder){
				bwDom1.write("\t");
				bwDom1.write(dom);
			}
			bwDom1.write("\n");
			for(String gene : id2net.get(id).keySet()){
				bwDom1.write(gene);
				for(String dom : domainOrder){
					bwDom1.write("\t");
					if(gene2domains1.get(gene).contains(dom)){
						bwDom1.write("1");
					}
					else if(gene2domains2.get(gene).contains(dom)){
						bwDom1.write("2");
					}
					else{
						bwDom1.write("0");
					}
				}
				bwDom1.write("\n");
			}
			bwDom1.flush();
			bwDom1.close();
			
			if(selected.contains(id)){
				BufferedWriter bwDom = new BufferedWriter(new FileWriter(outDir+"domains\\"+id+"_domains_annotated.tab"));
				bwDom.write("gene");
				for(String dom : domainOrder){
					bwDom.write("\t");
					bwDom.write(dom);
				}
				bwDom.write("\n");
				for(String gene : id2net.get(id).keySet()){
					if(!id2inphi.get(id).contains(gene)){
						continue;
					}
					bwDom.write(gene);
					for(String dom : domainOrder){
						bwDom.write("\t");
						if(gene2domains1.get(gene).contains(dom)){
							bwDom.write("1");
						}
						else if(gene2domains2.get(gene).contains(dom)){
							bwDom.write("2");
						}
						else{
							bwDom.write("0");
						}
					}
					bwDom.write("\n");
				}
				bwDom.flush();
				bwDom.close();
			}
		}
		
		
		
		//id2inphi
		
		//Domains feature
		//domainCount.add(domain);
		//for(String id : id2net.get(sel).get(gene)){}


	}
}
