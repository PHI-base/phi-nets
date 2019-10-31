package uk.ac.rothamsted.phinets;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.collect.SetMultimap;

import uk.ac.rothamsted.phinets.util.DefaultSuppliers;
import uk.ac.rothamsted.phinets.util.DelimitedIO;

public class DataCleaningScripts {

	public static void main(String[] args) throws Exception {
		//buildNetwork();
		
		//stats();
		//buildDdi();
		//processKbdock();
		//domineDdi();
		//buildDdi();
		//domineDdi();
		//clean();
		//parse3did();
		
		//condenceOma();
		//interologNetworkOma();
		//interologNetworkOrthoDb();
		//interologNetworkOrthoDb();
		
		/*
		clean();
		buildDdi();
		String dir = "";
		String file3did = "";
		*/
		
		//processBlast();
		 overlapStats();
		
	}
	
	public static void overlapStats() throws Exception{
		String dir = "D:\\projects\\path_net_combined\\networks\\complete\\";
		Map<String, DescriptiveStatistics> counts = new HashMap<String, DescriptiveStatistics>();
		Map<String, DescriptiveStatistics> percent = new HashMap<String, DescriptiveStatistics>();
		//Set<String> ids = new HashSet<String>(Arrays.asList(new String[]{"","","","","",""}));
		for(File f : new File(dir).listFiles()){
			BufferedReader br = new BufferedReader(new FileReader(f));
			String line = br.readLine();
			String name = f.getName().split("_interactome")[0];
			HashMultiset<String> sm = HashMultiset.create();
			double total = 0d;
			while ((line = br.readLine()) != null) {
				Set<String> tmp = new HashSet<String>();
				String[] data = line.split("\t");
				for(String d : data[4].split(",")){
					tmp.add(d);
				}
				for(String d : data[5].split(",")){
					tmp.add(d);
				}
				tmp.remove("");
				if(!data[4].trim().isEmpty()){
					tmp.add("ddi");
				}
				if(!data[5].trim().isEmpty()){
					tmp.add("interolog");
				}
				//sm.add(genId(tmp));
				sm.addAll(tmp);
				total+=1d;
			}
			//System.err.print(name);
			for(com.google.common.collect.Multiset.Entry<String> ent : sm.entrySet()){
				DescriptiveStatistics ds = counts.get(ent.getElement());
				if(ds == null){
					ds = new  DescriptiveStatistics();
					counts.put(ent.getElement(), ds);
				}
				//ds.addValue((((double)ent.getCount()) /total)*100d);
				ds.addValue((double)ent.getCount());
				
				ds = percent.get(ent.getElement());
				if(ds == null){
					ds = new  DescriptiveStatistics();
					percent.put(ent.getElement(), ds);
				}
				//ds.addValue((((double)ent.getCount()) /total)*100d);
				ds.addValue(((double)ent.getCount())/total);
			}
			//System.err.print("\n");
			br.close();
		}
		for(Entry<String, DescriptiveStatistics> ent: counts.entrySet()){
			//double stderr = ent.getValue().getStandardDeviation()/Math.sqrt(ent.getValue().getN());
			//System.err.println(ent.getKey()+"\t"+ent.getValue().getMean()+"\t"+stderr);
			System.err.println(ent.getKey()+"\t"+ent.getValue().getSum()+"\t"+percent.get(ent.getKey()).getMin()+"\t"+percent.get(ent.getKey()).getMax());
		}
	}
	
	public static String genId(Set<String> set){
		List<String> list = new ArrayList<String>(set);
		Collections.sort(list);
		StringBuffer sb = new StringBuffer();
		for(String s :list){
			if(sb.length() != 0){
				sb.append("&");
			}
			sb.append(s);
		}
		return sb.toString();
	}
	
	
	public static void processBlast() throws Exception{
		String inFile = "D:\\projects\\path_net\\seqsim\\decypher_file\\blastP_allFungiVSallFungi_out_050516_final_01.txt";
		String dir = "D:\\projects\\path_net\\seqsim\\fa\\";
		String []  exclude = new String[]{"Aspergillus_fumigatus_z5.ASM102932v1.31.pep.all.fa","Verticillium_dahliaejr2.GCA_000400815.2.31.pep.all.fa","Aspergillus_fumigatusa1163.CADRE.31.pep.all.fa"};
		String outFile = "D:\\projects\\path_net\\seqsim\\decypher_file\\blastP_allFungiVSallFungi_out_050516_final_01[removed_strains].txt";
		
		Set<String> remove = new HashSet<String>();
		for(String ex : exclude){
			BufferedReader br = new BufferedReader(new FileReader(dir+ex));
			String line;
			while ((line = br.readLine()) != null) {
				if(!line.startsWith(">")){
					continue;
				}
				String[] data = line.split(" ");
				String v = data[0].replace(">", "").trim();
				remove.add(v);
				System.err.println(v);
			}
			br.close();
		}
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		BufferedReader br = new BufferedReader(new FileReader(inFile));
		String line;
		int counter = 0 ;
		while ((line = br.readLine()) != null) {
			String [] data = line.split("\t");
			String from = data[1].split("\\|")[0].trim();
			String to = data[4].split("\\|")[0].trim();
			if(remove.contains(from) || remove.contains(to)){
				counter++;
				continue;
			}
			bw.write(line+"\n");
			
		}
		br.close();
		bw.flush();
		bw.close();
		
	}
	
	public static void parse3did() throws Exception{
		String file = "D:\\projects\\path_net\\3did_flat\\3did_flat_Jun_29_2015.dat";
		String out = "D:\\projects\\path_net\\3did_flat\\3did.tab";
		Pattern p = Pattern.compile("(PF[0-9]+).*?(PF[0-9]+)");
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		while ((line = br.readLine()) != null) {
			if(line.startsWith("#=ID")){
				Matcher m = p.matcher(line);
				if(m.find()){
					String a = m.group(1);
					String b = m.group(2);
					bw.write(a+"\t"+b+"\n");
				}
				else{
					System.err.println(line);
				}
			}
		}
		br.close();
		bw.flush();
		bw.close();
	}
	
	public static void interologNetworkOrthoDb() throws Exception{
		String sw = "D:\\projects\\path_net\\orthologs\\fg2swiss.txt";
		String tr = "D:\\projects\\path_net\\orthologs\\fg2trembl.txt";
		SetMultimap<String, String> uni2fg = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		DelimitedIO.readSimpleMap(uni2fg, sw, "\t", 1, 0);
		DelimitedIO.readSimpleMap(uni2fg, tr, "\t", 1, 0);
		Set<String> fgUniprotIds = new HashSet<String>(uni2fg.keySet());
		
		String orthodb = "D:\\projects\\path_net\\orthoDb\\ODB8_EukOGs_genes_ALL_levels.txt";
		BufferedReader br = new BufferedReader(new FileReader(orthodb));
		String line = br.readLine();
		SetMultimap<String, String> u2o = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		SetMultimap<String, String> o2u = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		Set<String> os = new HashSet<String>();
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			String u = data[2].split("_")[0].trim();
			if(fgUniprotIds.contains(u)){
				os.add(data[1].trim());
			}
		}
		br.close();
		System.err.println("Pass one finished.");
		System.err.println("Orthologous groups: "+os.size());
		
		br = new BufferedReader(new FileReader(orthodb));
		br.readLine();
		Set<String> set1 = new HashSet<String>();
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			String u = data[2].split("_")[0].trim();
			String o = data[1].trim();
			if(os.contains(o)){
				u2o.put(u, o);
				set1.add(u);
				if(fgUniprotIds.contains(u)){
					o2u.put(o, u);					
				}
			}
		}
		br.close();
		System.err.println("Pass two finished.");
		System.err.println(o2u.keySet().size());
		System.err.println(u2o.keySet().size());
		
		Set<String> fungiTaxids = new HashSet<String>(){{
			addAll(DelimitedIO.readColumn("D:\\projects\\path_net\\intact\\taxonomy_all_fungi.txt", 0, "\t"));	
		}};

		String intact = "D:\\projects\\path_net\\intact\\intact.txt";
		//String file = "D:\\projects\path_net\\intact\\biogrid_head.tab";
		br = new BufferedReader(new FileReader(intact));
		line = br.readLine();
		SetMultimap<String, String> interaction2taxa = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		SetMultimap<String, String> taxa2interaction = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		Map<String, Double> edges = new HashMap<String, Double>();
		SetMultimap<String, Integer> times = Multimaps.newSetMultimap(Maps.<String, Collection<Integer>> newHashMap(), DefaultSuppliers.<Integer> set());
		int i = 0;
		int j = 0;
		Set<String> iu = new HashSet<String>();
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			try{
				String t1 = data[9].split("\\:")[1].split("\\(")[0];
				String t2 = data[10].split("\\:")[1].split("\\(")[0];
				if(!t1.equals(t2)){
					continue;
				}
				if(data[35].equals("TRUE")){
					continue;
				}
				if(!data[0].contains("uniprotkb:") || !data[1].contains("uniprotkb:")){
					continue;
				}

				String ua = data[0].split("\\:")[1].split("-")[0];
				Set<String> aFg = new HashSet<String>();
				String ub = data[1].split("\\:")[1].split("-")[0];
				Set<String> bFg = new HashSet<String>();
				if(ua.equals(ub)){
					continue;
				}
				
				/*if(t1.equals("5518")){
					for(String a : uni2fg.get(ua)){
						for(String b : uni2fg.get(ub)){
							String pair;
							if(a.equals(b)){
								continue;
							}
							if(a.compareTo(b) > 0){
								pair = a +"\t" +b;
							}
							else{
								pair = b +"\t" +a;
							}
							System.err.println(pair);
							interaction2taxa.put(pair, t1);
							taxa2interaction.put(t1, pair);
							times.put(pair, i);
							Double conf = edges.get(pair);
							Double conf2 = Double.NaN;
							try{
								conf2 = Double.valueOf(data[14].split("\\:")[1]);
							}
							catch(java.lang.NumberFormatException e){
								conf2 = Double.NaN;
							}

							if(conf == null || conf2 > conf){
								edges.put(pair, conf2);
							}
						}
					}
					continue;
				}
				*/
				
				//if(!u2o.get(ua).isEmpty() && !u2o.get(ub).isEmpty() ){
				//	i++;
				//}
				//else{
					//System.err.println(ua + " - "+ub);
				//}
				//
				iu.add(ua);
				iu.add(ub);
				
				for(String o : u2o.get(ua)){
					for(String uni : o2u.get(o)){
						aFg.addAll(uni2fg.get(uni));
					}
				}
				
				for(String o : u2o.get(ub)){
					for(String uni : o2u.get(o)){
						bFg.addAll(uni2fg.get(uni));
					}
				}
				
				for(String a : aFg){
					for(String b : bFg){
						String pair;
						if(a.equals(b)){
							continue;
						}
						if(a.compareTo(b) > 0){
							pair = a +"\t" +b;
						}
						else{
							pair = b +"\t" +a;
						}
						times.put(pair, i);
						Double conf = edges.get(pair);
						Double conf2 = Double.NaN;
						try{
							conf2 = Double.valueOf(data[14].split("\\:")[1]);
						}
						catch(java.lang.NumberFormatException e){
							conf2 = Double.NaN;
						}

						if(conf == null || conf2 > conf){
							edges.put(pair, conf2);
						}
						
						interaction2taxa.put(pair, t1);
						taxa2interaction.put(t1, pair);
					}
				}
				
			}
			catch(java.lang.ArrayIndexOutOfBoundsException e){
				continue;
			}
			//i++;
		}
		br.close();
		System.err.println("All: "+iu.size());
		iu.retainAll(set1);
		System.err.println("Common: "+iu.size());
		
		System.err.println("Writing file");
		
		BufferedWriter bw = new BufferedWriter(new FileWriter("D:\\projects\\path_net\\fg_interactome[orthoDb_all].tab"));
		for(Entry<String, Collection<String>> ent : interaction2taxa.asMap().entrySet()){
			String flag = "0";
			if(ent.getValue().contains("5518")){
				flag = "2";
			}
			else {
				for(String tid : ent.getValue()){
					if(fungiTaxids.contains(tid)){
						flag = "1";
						break;
					}
				}
			}
			
			bw.write(ent.getKey()+"\t"+ent.getValue().size()+"\t"+flag+"\t"+String.valueOf(edges.get(ent.getKey()))+"\t"+String.valueOf(times.get(ent.getKey()).size())+"\n");
		}
		bw.flush();
		bw.close();
		
		System.err.println("Finished");
	}

	public static void interologNetworkOma() throws Exception{
		String sw = "D:\\projects\\path_net\\orthologs\\fg2swiss.txt";
		String tr = "D:\\projects\\path_net\\orthologs\\fg2trembl.txt";
		String oma = "D:\\projects\\path_net\\oma\\fg_oma.tab";
		String omau = "D:\\projects\\path_net\\oma\\oma-uniprot.txt";
		SetMultimap<String, String> uni2fg = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		DelimitedIO.readSimpleMap(uni2fg, sw, "\t", 1, 0);
		DelimitedIO.readSimpleMap(uni2fg, tr, "\t", 1, 0);
		Set<String> fgUniprotIds = new HashSet<String>(uni2fg.keySet());
		Set<String> fgOmaIds = new HashSet<String>();
		Set<String> allOmaIds = new HashSet<String>();
		allOmaIds.addAll(DelimitedIO.readColumn(oma, 0, "\t"));
		allOmaIds.addAll(DelimitedIO.readColumn(oma, 1, "\t"));
		SetMultimap<String, String> oma2uni = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		SetMultimap<String, String> uni2oma = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		SetMultimap<String, String> oma2fgOma = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		BufferedReader br = new BufferedReader(new FileReader(omau));
		String line;
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			if(allOmaIds.contains(data[0])){
				oma2uni.put(data[0], data[1]);
				uni2oma.put(data[1], data[0]);
				if(fgUniprotIds.contains(data[1])){
					fgOmaIds.add(data[0]);
				}
			}
			
		}
		br.close();
		System.err.println(fgUniprotIds.size());
		Set<String> uniprot = new HashSet<String>(DelimitedIO.readColumn(omau, 1, "\t"));
		System.err.println(uniprot.size());
		uniprot.retainAll(fgUniprotIds);
		System.err.println(uniprot.size());
		
		System.err.println("Finished reading uniprot mappings");
		br = new BufferedReader(new FileReader(oma));
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			if(fgOmaIds.contains(data[0])){
				oma2fgOma.put(data[1], data[0]);
			}
			else{
				oma2fgOma.put(data[0], data[1]);
			}
		}
		br.close();
		System.err.println("Finished parsing oma orthology mappings");

		Set<String> fungiTaxids = new HashSet<String>(){{
			addAll(DelimitedIO.readColumn("D:\\projects\\path_net\\intact\\taxonomy_all_fungi.txt", 0, "\t"));	
		}};
		
		String intact = "D:\\projects\\path_net\\intact\\intact.txt";
		//String file = "D:\\projects\path_net\\intact\\biogrid_head.tab";
		br = new BufferedReader(new FileReader(intact));
		line = br.readLine();
		SetMultimap<String, String> interaction2taxa = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		SetMultimap<String, String> taxa2interaction = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		Map<String, Double> edges = new HashMap<String, Double>();
		SetMultimap<String, Integer> times = Multimaps.newSetMultimap(Maps.<String, Collection<Integer>> newHashMap(), DefaultSuppliers.<Integer> set());
		int i = 0;
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			try{
				String t1 = data[9].split("\\:")[1].split("\\(")[0];
				String t2 = data[10].split("\\:")[1].split("\\(")[0];
				if(!t1.equals(t2)){
					continue;
				}
				if(data[35].equals("TRUE")){
					continue;
				}
				if(!data[0].contains("uniprotkb:") || !data[1].contains("uniprotkb:")){
					continue;
				}

				String ua = data[0].split("\\:")[1];
				Set<String> aFg = new HashSet<String>();
				String ub = data[1].split("\\:")[1];
				Set<String> bFg = new HashSet<String>();
				if(ua.equals(ub)){
					continue;
				}

				if(t1.equals("5518")){
					for(String a : uni2fg.get(ua)){
						for(String b : uni2fg.get(ub)){
							String pair;
							if(a.compareTo(b) > 0){
								pair = a +"\t" +b;
							}
							else{
								pair = b +"\t" +a;
							}
							interaction2taxa.put(pair, t1);
							taxa2interaction.put(t1, pair);
							times.put(pair, i);
							Double conf = edges.get(pair);
							Double conf2 = Double.NaN;
							try{
								conf2 = Double.valueOf(data[14].split("\\:")[1]);
							}
							catch(java.lang.NumberFormatException e){
								conf2 = Double.NaN;
							}

							if(conf == null || conf2 > conf){
								edges.put(pair, conf2);
							}
						}
					}
					continue;
				}
				
				
				for(String o : uni2oma.get(ua)){
					for(String fgo : oma2fgOma.get(o)){
						for(String uni : oma2uni.get(fgo)){
							aFg.addAll(uni2fg.get(uni));
						}
					}
				}
				
				for(String o : uni2oma.get(ub)){
					for(String fgo : oma2fgOma.get(o)){
						for(String uni : oma2uni.get(fgo)){
							bFg.addAll(uni2fg.get(uni));
						}
					}
				}
				
				for(String a : aFg){
					for(String b : bFg){
						String pair;
						if(a.compareTo(b) > 0){
							pair = a +"\t" +b;
						}
						else{
							pair = b +"\t" +a;
						}
						times.put(pair, i);
						Double conf = edges.get(pair);
						Double conf2 = Double.NaN;
						try{
							conf2 = Double.valueOf(data[14].split("\\:")[1]);
						}
						catch(java.lang.NumberFormatException e){
							conf2 = Double.NaN;
						}

						if(conf == null || conf2 > conf){
							edges.put(pair, conf2);
						}
						
						interaction2taxa.put(pair, t1);
						taxa2interaction.put(t1, pair);
					}
				}
				
			}
			catch(java.lang.ArrayIndexOutOfBoundsException e){
				continue;
			}
			i++;
		}
		br.close();
		
		System.err.println("Writing file");
		
		BufferedWriter bw = new BufferedWriter(new FileWriter("D:\\projects\\path_net\\fg_interactome[oma_all].tab"));
		for(Entry<String, Collection<String>> ent : interaction2taxa.asMap().entrySet()){
			String flag = "0";
			if(ent.getValue().contains("5518")){
				flag = "2";
			}
			else {
				for(String tid : ent.getValue()){
					if(fungiTaxids.contains(tid)){
						flag = "1";
						break;
					}
				}
			}
			
			bw.write(ent.getKey()+"\t"+ent.getValue().size()+"\t"+flag+"\t"+String.valueOf(edges.get(ent.getKey()))+"\t"+String.valueOf(times.get(ent.getKey()).size())+"\n");
		}
		bw.flush();
		bw.close();
		
		System.err.println("Finished");
	}

	public static void condenceOma() throws Exception{
		String sw = "D:\\projects\\path_net\\orthologs\\fg2swiss.txt";
		String tr = "D:\\projects\\path_net\\orthologs\\fg2trembl.txt";
		String oma = "D:\\projects\\path_net\\oma\\oma-pairs.txt";
		String omau = "D:\\projects\\path_net\\oma\\oma-uniprot.txt";
		SetMultimap<String, String> uni2fg = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		DelimitedIO.readSimpleMap(uni2fg, sw, "\t", 1, 0);
		DelimitedIO.readSimpleMap(uni2fg, tr, "\t", 1, 0);
		SetMultimap<String, String> oma2uni = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		SetMultimap<String, String> oma2uni2 = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		BufferedReader br = new BufferedReader(new FileReader(omau));
		String line;
		while ((line = br.readLine()) != null) {
			if(line.startsWith("#")){
				continue;
			}
			String[] data = line.split("\t");
			if(uni2fg.containsKey(data[1])){
				oma2uni.put(data[0], data[1]);
			}
		}
		br.close();
		System.err.println("Parsed uniprot mappings.");		
		
		Set<String> set = oma2uni.keySet();
		br = new BufferedReader(new FileReader(oma));
		
		BufferedWriter bw = new BufferedWriter(new FileWriter("D:\\projects\\path_net\\fg_oma.tab"));
		while ((line = br.readLine()) != null) {
			if(line.startsWith("#")){
				continue;
			}
			String[] data = line.split("\t");
			if(set.contains(data[0]) || set.contains(data[1])){
				bw.write(line+"\n");
			}
		}
		bw.flush();
		bw.close();
	}
	
	public static void clean() throws Exception{
		//DelimitedIO.removeBlanks("D:\\projects\\path_net\\go\\mart_export (1).txt", "D:\\projects\\path_net\\go\\go_fgsg.tab", true);
		//DelimitedIO.removeBlanks("D:\\projects\\path_net\\go\\mart_export.txt", "D:\\projects\\path_net\\go\\go_trid.tab", true);
		//DelimitedIO.removeBlanks("D:\\projects\\path_net\\phi-base_biomart\\mart_export (2).txt", "D:\\projects\\path_net\\phi-base_biomart\\trid_phibase.tab", true);
		DelimitedIO.removeBlanks("D:\\projects\\path_net\\orthologs\\mart_export.txt", "D:\\projects\\path_net\\orthologs\\fg2swiss.txt", true);
		DelimitedIO.removeBlanks("D:\\projects\\path_net\\orthologs\\mart_export (1).txt", "D:\\projects\\path_net\\orthologs\\fg2trembl.txt", true);
	}

	public static void domineDdi() throws Exception{
		String dFile = "D:\\projects\\path_net\\domine_2010\\INTERACTION.txt";
		BufferedReader br = new BufferedReader(new FileReader(dFile));
		String line;
		Set<String> set = new HashSet<String>();
		BufferedWriter bw = new BufferedWriter(new FileWriter("D:\\projects\\path_net\\domine_2010\\domine_hc.tab"));

		while ((line = br.readLine()) != null) {
			String[] data = line.split("\\|");
			String pair = data[0]+"\t"+data[1]+"\n";
			if(data[1].compareTo(data[0]) > 0){
				pair = data[1]+"\t"+data[0]+"\n";
			}
			if(data[17].equals("HC")){ // || data[17].equals("LC")){
				if(!set.contains(pair)){
					set.add(pair);
					bw.write(pair);
				}
			}

		}
		br.close();
		bw.flush();
		bw.close();
	}
	
	public static void processKbdock()throws Exception{
		Pattern p = Pattern.compile("(PF[0-9]*)_.*?(PF[0-9]*)_");
		String dir = "D:\\projects\\path_net\\KBDOCK\\zip_sup_rep_ddi_intra_homo\\";
		BufferedWriter bw = new BufferedWriter(new FileWriter("D:\\projects\\path_net\\kbdocks_homo.tab"));
		
		Set<String> set = new HashSet<String>();
		for(File sub : new File(dir).listFiles()){
			if(sub.isDirectory()){
				for(File f : sub.listFiles()){
					Matcher m = p.matcher(f.getName());
					if(m.find()){
						String line = m.group(1).trim()+"\t"+m.group(2).trim();
						if(!set.contains(line)){
							set.add(line);
							bw.write(line+"\n");
						}

					}
					continue;
				}
			}
		}
		bw.flush();
		bw.close();
			
	}
	
	public static void buildDdi() throws Exception{
		String dFile = "D:\\projects\\path_net\\anno\\output-correctedAfterManualOverlappingSolving-GOLDfile.txt";
		
		/*
		String ddFile = "D:\\projects\\path_net\\iPFAM\\heterodomain_interaction.csv";
		List<String[]> pairs = new ArrayList<String[]>();
		BufferedReader br = new BufferedReader(new FileReader(ddFile));
		String line = br.readLine();
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			pairs.add(new String[]{data[0], data[2]});
		}
		br.close();
		*/
		
		//String ddFile = "D:\\projects\\path_net\\domine_2010\\domine_all.tab";
		//String ddFile = "D:\\projects\\path_net\\kbdocks_all.tab";
		String ddFile = "D:\\projects\\path_net\\3did_flat\\3did.tab";
		List<String[]> pairs = new ArrayList<String[]>();
		BufferedReader br = new BufferedReader(new FileReader(ddFile));
		String line;// = br.readLine();
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			pairs.add(new String[]{data[0], data[1]});
			
		}
		br.close();
		
		br = new BufferedReader(new FileReader(dFile));
		SetMultimap<String, String> map = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			String dom = data[1].split("\\.")[0];
			map.put(dom, data[0].trim());
		}
		br.close();

		HashMultiset<String> ms = HashMultiset.create();
		for(String [] pair : pairs){
			//System.err.println(pair[0]+"\t"+pair[1]);
			Collection<String> as = map.get(pair[0]);
			Collection<String> bs = map.get(pair[1]);
			for(String a :as){
				for(String b : bs){
					if(!a.equals(b)){
						if(a.compareTo(b) > 0){
							ms.add(a+"\t"+b);
						}
						else{
							ms.add(b+"\t"+a);
						}
					}
				}
			}
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter("D:\\projects\\path_net\\fg_ddi_[3did].tab"));
		for(com.google.common.collect.Multiset.Entry<String> en : ms.entrySet()){
			bw.write(en.getElement()+"\t"+String.valueOf(en.getCount()+"\n"));
		}
		bw.flush();
		bw.close();
	}
	
	public static void stats()throws Exception{
		Set<String> ids = new HashSet<String>(){{
			addAll(DelimitedIO.readColumn("D:\\projects\\path_net\\intact\\taxonomy_all_fungi.txt", 0, "\t"));	
		}};
		
		//String file = "D:\\projects\\path_net\\intact\\fungi.txt";
		String file = "D:\\projects\\path_net\\intact\\intact.txt";
		//String file = "D:\\projects\path_net\\intact\\biogrid_head.tab";
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		HashMultiset<String> times = HashMultiset.create();
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			try{
				String t1 = data[9].split("\\:")[1].split("\\(")[0];
				String t2 = data[10].split("\\:")[1].split("\\(")[0];
				if(t1.equals(t2) && ids.contains(t1) && ids.contains(t2)){
					times.add(t1);
				}	
			}
			catch(java.lang.ArrayIndexOutOfBoundsException e){
				continue;
			}

		}
		br.close();
		BufferedWriter bw = new BufferedWriter(new FileWriter("D:\\projects\\path_net\\intact\\fungi_stats.txt"));
		for(com.google.common.collect.Multiset.Entry<String> ent : times.entrySet()){
			//System.err.println(ent.getElement()+"\t"+ent.getCount());
			bw.write(ent.getElement()+"\t"+ent.getCount()+"\n");
		}
		bw.flush();
		bw.close();
		System.err.println("done");
	}
	
	public static final void buildNetwork()throws Exception{
		String dir = "D:\\projects\\path_net\\";
		/* sp
		String species = "sp";
		String orthFile = dir+"orthologs\\fg2sp.txt";
		String uni1File = dir+"orthologs\\sp_swiss.txt";
		String uni2File = dir+"orthologs\\sp_trembl.txt";
		String intactFile = dir+"intact\\Schizosacc.txt";
		 */
		
		/* sc
		String species = "sc";
		String orthFile = dir+"orthologs\\fg2sc.txt";
		String uni1File = dir+"orthologs\\sc_swiss.txt";
		String uni2File = dir+"orthologs\\sc_termbl.txt";
		String intactFile = dir+"intact\\Saccharomy.txt";
		*/
		
		/* af */
		String species = "af";
		String orthFile = dir+"orthologs\\fg2af.txt";
		String uni1File = dir+"orthologs\\af_swiss.txt";
		String uni2File = dir+"orthologs\\af_trembl.txt";
		String intactFile = dir+"intact\\Aspergillu.txt";
		
		
		SetMultimap<String, String> om = DelimitedIO.readEnsemblMap(orthFile, 1, 0, null, 2);
		System.err.println("Orthologs: "+om.size());
		SetMultimap<String, String> um = DelimitedIO.readEnsemblMap(uni1File, 1, 0, null, 2);
		SetMultimap<String, String> u2m = DelimitedIO.readEnsemblMap(uni2File, 1, 0, null, 2);
		um.putAll(u2m);
		System.err.println("Uniprot mappings: "+um.size());
		
		
		BufferedReader br = new BufferedReader(new FileReader(intactFile));
		String line = br.readLine();
		Set<String> nodes = new HashSet<String>();
		Map<String,Double> edges = new HashMap<String, Double>();
		HashMultiset<String> times = HashMultiset.create();
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			if(data[35].equals("TRUE")){
				continue;
			}
			if(!data[0].contains("uniprotkb:") || !data[1].contains("uniprotkb:")){
				continue;
			}
			String ua = data[0].split("\\:")[1];
			List<String> as = getAll(um.get(ua), om);
			String ub = data[1].split("\\:")[1];
						
			List<String> bs = getAll(um.get(ub), om);
			for(String a : as){
				for(String b : bs){
					
					if(a.equals(b)){
						continue;
					}
					String key;
					if(a.compareTo(b) < 0){
						key = a +"\t"+ b;
					}
					else{
						key =  b +"\t"+ a;
					}
					
					Double conf = edges.get(key);
					Double conf2 = Double.NaN;
					try{
						conf2 = Double.valueOf(data[14].split("\\:")[1]);
					}
					catch(java.lang.NumberFormatException e){
						conf2 = Double.NaN;
						//System.err.println(line);
						//continue;
					}
					times.add(key);
					
					
					if(conf == null || conf2 > conf){
						
						edges.put(key, conf2);
						nodes.add(a);
						nodes.add(b);
					}


				}
			}
		}
		br.close();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(dir+"fg_interactome["+species+"].tab"));
		for(Entry<String, Double> ent: edges.entrySet()){
			bw.write(ent.getKey()+"\t"+ent.getValue().toString()+"\t"+String.valueOf(times.count(ent.getKey()))+"\n");
		}
		bw.flush();
		bw.close();
		
	}
	
	public static final List<String> getAll(Collection<String> col, Multimap<String, String> map){
		Set<String> set = new HashSet<String>();
		for(String c : col){
			set.addAll(map.get(c));
		}
		return new ArrayList<String>(set);
	}
}
