package uk.ac.rothamsted.phinets;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sourceforge.OntologyGraph;
import net.sourceforge.OntologyGraph.Entity;
import net.sourceforge.OntologyGraph.Term;
import net.sourceforge.utils.ArrayKey;
import uk.ac.rothamsted.phinets.util.DefaultSuppliers;
import uk.ac.rothamsted.phinets.util.DelimitedIO;

import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.collect.SetMultimap;
import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;

public class InteractomeBuilder {
	
	public static List<String[]> readStageConfig(String file) throws Exception{
		List<String[]> list = new ArrayList<String[]>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			list.add(data);
		}
		br.close();
		return list;
	}
	
	public static void mergeGoStats() throws Exception{
		String baseDir = "D:\\projects\\path_net_combined\\goa_stats\\";
		String [] types = new String[]{"bp", "mf", "cc"};
		Map<String, Integer> subtypes = new HashMap<String, Integer>(){{
			put("kbdocks_all", 0);
			put("domine_hc", 1);
			put("3did", 2);
			put("sc", 3);
			put("sp", 4);
			put("rnd", 5);
		}};
		
		for(String type : types){
			BufferedWriter bw = new BufferedWriter(new FileWriter(baseDir+type+"_all.tab"));
			boolean first = true;
			for(File file: new File(baseDir+type).listFiles()){
				BufferedReader br = new BufferedReader(new FileReader(file));
				String line = br.readLine();
				String[] data = line.split("\t");
				if(first){
					bw.write(data[0]+"\t"+data[1]+"\t"+data[2]+"\t"+data[4]+"\tkbdocks_all\tdomine_hc\t3did\tsc\tsp\trnd\n");
					first = false;
				}
				
				while ((line = br.readLine()) != null) {
					String [] member = new String[]{"0","0","0","0","0","0"};
					data = line.split("\t");
					String [] tt = data[3].split(", ");
					for(String t : tt){
						member[subtypes.get(t)] = "1";
					}
					bw.write(data[0]+"\t"+data[1]+"\t"+data[2]+"\t"+data[4]);
					for(String m : member){
						bw.write("\t"+m);
					}
					bw.write("\n");
				}
				
				br.close();
			}
			for(File file: new File(baseDir+type+"_rnd").listFiles()){
				BufferedReader br = new BufferedReader(new FileReader(file));
				String line = br.readLine();
				while ((line = br.readLine()) != null) {
					String [] member = new String[]{"0","0","0","0","0","0"};
					String[] data = line.split("\t");
					String [] tt = data[3].split(", ");
					for(String t : tt){
						member[subtypes.get(t)] = "1";
					}
					bw.write(data[0]+"\t"+data[1]+"\t"+data[2]+"\t"+data[4]);
					for(String m : member){
						bw.write("\t"+m);
					}
					bw.write("\n");
				}
				
				br.close();
			}
			bw.flush();
			bw.close();
		}
	}
	
	public static void overlapStats(String dir) throws Exception{
		Set<String> set = new HashSet<String>(){{add("kbdocks_all");add("domine_hc");add("3did");add("sc");add("sp");}};
		Set<Set<String>> combos = Sets.powerSet(set);
		HashMultiset<Set<String>> counter = HashMultiset.create();
		Map<String, HashMultiset<Set<String>>> map = new HashMap<String, HashMultiset<Set<String>>>();
		for(File f : new File(dir).listFiles()){
			HashMultiset<Set<String>> subCounter = HashMultiset.create();
			String [] split = f.getName().split("_");
			map.put(split[0]+"_"+split[1], subCounter);
			BufferedReader br = new BufferedReader(new FileReader(f));
			String line = br.readLine();
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				Set<String> tmp = new HashSet<String>();
				data[4] = data[4].trim();
				data[5] = data[5].trim();
				if(!data[4].isEmpty()){
					tmp.addAll(Arrays.asList(data[4].split(",")));	
				}
				if(!data[5].trim().isEmpty()){
					tmp.addAll(Arrays.asList(data[5].split(",")));	
				}
				for(Set<String> combo : combos){
					if(tmp.size() == combo.size() && combo.containsAll(tmp)){
						counter.add(combo);
						subCounter.add(combo);
					}
				}
				
			}
			br.close();
		}
		
		for(com.google.common.collect.Multiset.Entry<Set<String>> ent: counter.entrySet()){
			System.err.println("Total\t"+makeLabel(ent.getElement())+"\t"+ent.getCount());
		}
		for(Entry<String, HashMultiset<Set<String>>> e1 : map.entrySet()){
			for(com.google.common.collect.Multiset.Entry<Set<String>> ent: e1.getValue().entrySet()){
				System.err.println(e1.getKey()+"\t"+makeLabel(ent.getElement())+"\t"+ent.getCount());
			}
		}
	}
	
	public static void goStats(String dir, String outputDir, String topDir) throws Exception{
		//String [] aspects = new String[]{"biological_process", "molecular_function"}; MicaEvaluator CompartmentEvaluator
		
		/*
		Map<String, double[]> stats = new HashMap<String, double[]>(){{
			put("kbdocks_all", new double[2]);
			put("domine_hc", new double[2]);
			put("3did", new double[2]);
			put("sc", new double[2]);
			put("sp", new double[2]);
		}};
		*/
		
		//String [] aspects = new String[]{"biological_process", "molecular_function"};
		
		for(File f : new File(dir).listFiles()){
			String [] split = f.getName().split("_");
			String species = split[0]+"_"+split[1];
			
			BufferedWriter bw_cmp = new BufferedWriter(new FileWriter(outputDir+"cc\\"+species+"_cc.tab"));
			BufferedWriter bw_cmp_rnd = new BufferedWriter(new FileWriter(outputDir+"cc_rnd\\"+species+"_cc_rnd.tab"));
			
			BufferedWriter bw_mf = new BufferedWriter(new FileWriter(outputDir+"mf\\"+species+"_mf.tab"));
			BufferedWriter bw_mf_rnd = new BufferedWriter(new FileWriter(outputDir+"mf_rnd\\"+species+"_mf_rnd.tab"));
			
			BufferedWriter bw_bp = new BufferedWriter(new FileWriter(outputDir+"bp\\"+species+"_bp.tab"));
			BufferedWriter bw_bp_rnd = new BufferedWriter(new FileWriter(outputDir+"bp_rnd\\"+species+"_bp_rnd.tab"));
			
			BufferedWriter [] bws = new BufferedWriter[]{bw_cmp, bw_cmp_rnd, bw_mf, bw_mf_rnd, bw_bp, bw_bp_rnd};

			for(BufferedWriter bw : bws){
				bw.write("id1\tid2\tspecies\tedge_sources\tvalue\n");
				bw.flush();
			}
			
			int []  counters = new int[4];
			
			MicaEvaluator me = new MicaEvaluator(topDir, species);
			CompartmentEvaluator ce = new CompartmentEvaluator(topDir, species);
			
			BufferedReader br = new BufferedReader(new FileReader(f));
			String line;
			Set<String> nodes = new HashSet<String>();
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				Set<String> tmp = new HashSet<String>();
				data[4] = data[4].trim();
				data[5] = data[5].trim();
				nodes.add(data[0]);
				nodes.add(data[1]);
				if(!data[4].isEmpty()){
					tmp.addAll(Arrays.asList(data[4].split(",")));	
				}
				if(!data[5].trim().isEmpty()){
					tmp.addAll(Arrays.asList(data[5].split(",")));	
				}
				counters[3]+=1;
				Double bpValue = me.getMicaScore("biological_process", data[0], data[1]);
				String lbl = makeLabel(tmp);
				if(bpValue != null){
					bw_bp.write(data[0]+"\t"+data[1]+"\t"+species+"\t"+lbl+"\t"+String.valueOf(bpValue)+"\n");
					counters[0]+=1;
				}
				Double mfValue = me.getMicaScore("molecular_function", data[0], data[1]);
				if(mfValue != null){
					bw_mf.write(data[0]+"\t"+data[1]+"\t"+species+"\t"+lbl+"\t"+String.valueOf(mfValue)+"\n");
					counters[1]+=1;
				}
				Integer ccValue = ce.eval(data[0], data[1]);
				if(ccValue != null){
					bw_cmp.write(data[0]+"\t"+data[1]+"\t"+species+"\t"+lbl+"\t"+String.valueOf(ccValue)+"\n");
					counters[2]+=1;
				}
			}
			br.close();
			System.err.println(species+"\t"+String.valueOf(counters[0])+"\t"+String.valueOf(counters[1])+"\t"+String.valueOf(counters[2])+"\t"+String.valueOf(counters[3]));
			
			Set<String> tmp = new HashSet<String>(nodes);
			tmp.retainAll(me.getAllAnnotated("biological_process"));
			List<String> list = new ArrayList<String>(tmp);
			Set<ArrayKey<String>> done = new HashSet<ArrayKey<String>>();
			Random rnd = new Random(System.currentTimeMillis());
			while(done.size() < counters[0]){
				String a = list.get(rnd.nextInt(list.size()));
				String b = null;
				do{
					b = list.get(rnd.nextInt(list.size()));
				}while(a.equals(b));
				ArrayKey<String> key = ArrayKey.createSortedKey(a, b);
				if(done.contains(key)){
					continue;
				}
				done.add(key);
				Double bpValue = me.getMicaScore("biological_process", a, b);
				bw_bp_rnd.write(a+"\t"+b+"\t"+species+"\trnd\t"+String.valueOf(bpValue)+"\n");
			}
			
			tmp = new HashSet<String>(nodes);
			tmp.retainAll(ce.getAnnotated());
			list = new ArrayList<String>(tmp);
			done = new HashSet<ArrayKey<String>>();
			rnd = new Random(System.currentTimeMillis());
			while(done.size() < counters[2]){
				String a = list.get(rnd.nextInt(list.size()));
				String b = null;
				do{
					b = list.get(rnd.nextInt(list.size()));
				}while(a.equals(b));
				ArrayKey<String> key = ArrayKey.createSortedKey(a, b);
				if(done.contains(key)){
					continue;
				}
				done.add(key);
				Integer bpValue = ce.eval(a, b);
				bw_cmp_rnd.write(a+"\t"+b+"\t"+species+"\trnd\t"+String.valueOf(bpValue)+"\n");
			}
			
			tmp = new HashSet<String>(nodes);
			tmp.retainAll(me.getAllAnnotated("molecular_function"));
			list = new ArrayList<String>(tmp);
			done = new HashSet<ArrayKey<String>>();
			rnd = new Random(System.currentTimeMillis());
			while(done.size() < counters[1]){
				String a = list.get(rnd.nextInt(list.size()));
				String b = null;
				do{
					b = list.get(rnd.nextInt(list.size()));
				}while(a.equals(b));
				ArrayKey<String> key = ArrayKey.createSortedKey(a, b);
				if(done.contains(key)){
					continue;
				}
				done.add(key);
				Double bpValue = me.getMicaScore("molecular_function", a, b);
				bw_mf_rnd.write(a+"\t"+b+"\t"+species+"\trnd\t"+String.valueOf(bpValue)+"\n");
			}
			
			
			for(BufferedWriter bw : bws){
				bw.flush();
				bw.close();
			}
		}
	
	}
	
	
	public static String makeLabel(Collection<String> set){
		StringBuffer sb = new StringBuffer();
		boolean first = true;
		for(String s : set){
			if(!first){
				sb.append(", ");
			}
			sb.append(s);
			first = false;
		}
		return sb.toString();
	}

	public static String parseSc(String outDir, String shortId, String swissFile, String tremblFile, String intactFile, String species_prefix) throws Exception{
		SetMultimap<String, String> um = DelimitedIO.readEnsemblMap(swissFile, 1, 0, null, 2);
		SetMultimap<String, String> u2m = DelimitedIO.readEnsemblMap(tremblFile, 1, 0, null, 2);
		um.putAll(u2m);
		
		BufferedReader br = new BufferedReader(new FileReader(intactFile));
		String line = br.readLine();
		Set<String> nodes = new HashSet<String>();
		Map<String,Double> edges = new HashMap<String, Double>();
		HashMultiset<String> times = HashMultiset.create();
		Pattern num = Pattern.compile("intact-miscore:(.*)");
		
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			if(data[35].equals("TRUE")){
				continue;
			}
			if(!data[0].contains("uniprotkb:") || !data[1].contains("uniprotkb:")){
				continue;
			}
			String ua = data[0].split("\\:")[1];
			List<String> as = new ArrayList<String>(um.get(ua));
			String ub = data[1].split("\\:")[1];
						
			List<String> bs = new ArrayList<String>(um.get(ub));
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
					Matcher m = num.matcher(data[14]);
					if(m.find()){
						conf2 = Double.valueOf(m.group(1));
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
		
		String fileName = outDir+species_prefix+"_interactome["+shortId+"].tab";
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		for(Entry<String, Double> ent: edges.entrySet()){
			bw.write(ent.getKey()+"\t"+ent.getValue().toString()+"\t"+String.valueOf(times.count(ent.getKey()))+"\n");
		}
		bw.flush();
		bw.close();
		return fileName;
	}
	

	//String speciesListFile = "D:\\projects\\path_net_combined\\domains\\species_list.tab";
	//String orthoDir = "D:\\projects\\path_net_combined\\orthologs\\";
	//String goDir = "D:\\projects\\path_net_combined\\goa\\";
	//String phiDir = "D:\\projects\\path_net_combined\\phibase\\";
	public static void main(String[] args)throws Exception {
		//overlapStats("D:\\projects\\path_net_combined\\networks\\complete\\");
		//goStats("D:\\projects\\path_net_combined\\networks\\complete\\", "D:\\projects\\path_net_combined\\goa_stats\\", "D:\\projects\\path_net_combined\\") ;
		//mergeGoStats();
		System.err.println("Starting...");
		//System.exit(0);
		
		final String baseDir = "D:\\projects\\path_net_combined\\";
		final String outDir = baseDir+"networks\\";
		final String speciesList = baseDir+"species_list.tab";

		final String goaDir = baseDir+"goa\\";
		final String phiDir = baseDir+"phibase\\";
		final String goObo = baseDir+"go.obo";
		
		BufferedReader br = new BufferedReader(new FileReader(speciesList));
		String line;

		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			String species_prefix = data[0];
			
			String domainsFile = data[1];
			String ensemblId = data[2];
			List<String> networksByType = new ArrayList<String>();
			List<String> networksByDb = new ArrayList<String>();
			
			//Domain-domain interaction network
			final String ddiDir = baseDir+"ddi\\";
			final String domainsDir = baseDir+"domains\\";
			List<String> netDdiFiles = new ArrayList<String>();
			SetMultimap<String, ArrayKey<String>> map3 = Multimaps.newSetMultimap(Maps.<String, Collection<ArrayKey<String>>> newHashMap(), DefaultSuppliers.<ArrayKey<String>> set());
			for(File f : new File(ddiDir).listFiles()){
				String ddiSuffix = f.getName().split("\\.")[0];
				String ddiFile = f.getAbsolutePath();
				System.err.println("Processing "+species_prefix+" "+ddiSuffix);
				netDdiFiles.add(buildDdiNetwork(outDir, domainsDir+domainsFile, ddiFile, ddiSuffix, species_prefix, map3));
			}
			String fileName2 = outDir+"\\"+species_prefix+"_ddi_full_[combined].tab";
			BufferedWriter bw = new BufferedWriter(new FileWriter(fileName2));
			for(Entry<String, Collection<ArrayKey<String>>> en : map3.asMap().entrySet()){
				bw.write(en.getKey()+"\t"+String.valueOf(en.getValue().toString()+"\n"));
			}
			bw.flush();
			bw.close();
			
			System.err.println("Merging files...");
			String allDdiFile = mergeDdi(outDir, netDdiFiles, species_prefix);
			networksByType.add(allDdiFile);
			networksByDb.addAll(netDdiFiles);
			System.err.println("Done.");
			
			//Interolog networks construction for each interaction source species
			System.err.println("Building the interolog networks...");
			List<String> netInterologFiles = new ArrayList<String>();
			//orthoDir
			final String orthoDir = baseDir+"orthologs\\";
			final String idMappingsDir = baseDir+"mapping\\";
			final String intactDir = baseDir+"intact\\";
			for(File f : new File(orthoDir).listFiles()){
				if(!f.getName().startsWith(species_prefix)){
					continue;
				}
				String orthoFile = f.getAbsolutePath();
				String speciesSource = f.getName().split("_vs_")[1].split("\\.")[0];
				String shortId =  speciesSource.substring(0, 2);
				String swissFile = null;
				String tremblFile = null;
				for(File f1 : new File(idMappingsDir).listFiles()){
					//System.err.println(idMappingsDir+" ::" +f1.getName()+" :: "+shortId);
					if(!f1.getName().startsWith(shortId)){
						continue;
					}
					if(f1.getName().contains("swiss")){
						swissFile = f1.getAbsolutePath();
					}
					else if(f1.getName().contains("trembl")){
						tremblFile = f1.getAbsolutePath();
					}
				}
				String intactFile = null;
				for(File f1 : new File(intactDir).listFiles()){
					if(!f1.getName().startsWith(shortId)){
						continue;
					}
					intactFile = f1.getAbsolutePath();
					break;
				}
				
				System.err.println("Processing "+species_prefix+" "+shortId);
				if(f.getName().startsWith("Saccharomyces_cerevisiae_vs_scerevisiae")){
					netInterologFiles.add(parseSc(outDir, shortId, swissFile, tremblFile, intactFile, species_prefix));
					continue;
				}
				netInterologFiles.add(buildInterologNetwork(outDir, shortId, orthoFile, swissFile, tremblFile, intactFile, species_prefix));
			}
			System.err.println("Merging files...");
			String allInterologFile = mergeInterolog(outDir, netInterologFiles, species_prefix);
			networksByType.add(allInterologFile);
			networksByDb.addAll(netInterologFiles);
			System.err.println("Done.");
			String complete = buildCompleteNetwork(outDir+"complete\\", allInterologFile, allDdiFile, species_prefix);
			System.err.println("Complete network saved to :"+complete);
		}
		br.close();
	}
	
	private static final List<String[]> ddiConfig = new ArrayList<String[]>(){{
		add(new String[]{"output-correctedAfterManualOverlappingSolving-GOLDfile.txt", "ddi\\domine_hc.tab", "domine"});
		add(new String[]{"output-correctedAfterManualOverlappingSolving-GOLDfile.txt", "ddi\\kbdocks_all.tab", "kbdocks"});
		add(new String[]{"output-correctedAfterManualOverlappingSolving-GOLDfile.txt", "ddi\\3did.tab", "3did"});
	}};

	private static final List<String[]> interologConfig = new ArrayList<String[]>(){{
		add(new String[]{"orth_sp", "biomart\\fg2sp.txt", "biomart\\sp_swiss.txt", "biomart\\sp_trembl.txt", "intact\\sp.txt"});
		add(new String[]{"orth_sc", "biomart\\fg2sc.txt", "biomart\\sc_swiss.txt", "biomart\\sc_trembl.txt", "intact\\sc.txt"});
	}};
	
	
	
	public static void test(String[] args)throws Exception {
		final String workingDir = "D:\\projects\\path_net_combined\\";

		String species_prefix = "fg";
		
		//List<String> list = new ArrayList<String>(){{add(workingDir+"fg_interactome[ddi_all].tab");}};
		//goInteractomeValidation(workingDir, list, species_prefix);
		//System.exit(0);
		
		List<String> networksByType = new ArrayList<String>();
		List<String> networksByDb = new ArrayList<String>();
		
		System.err.println("Building the ddi networks...");
		//Domain-domain interaction networks construction for each domain interaction database
		List<String> netDdiFiles = new ArrayList<String>();
		for(String[] conf: ddiConfig){
			System.err.println("Processing "+species_prefix+" "+conf[2]);
			netDdiFiles.add(buildDdiNetwork(workingDir, workingDir+conf[0], workingDir+conf[1], conf[2], species_prefix, null));
		}
		System.err.println("Merging files...");
		String allDdiFile = mergeDdi(workingDir, netDdiFiles, species_prefix);
		networksByType.add(allDdiFile);
		networksByDb.addAll(netDdiFiles);
		System.err.println("Done.");
		
		//Interolog networks construction for each interaction source species
		System.err.println("Building the interolog networks...");
		List<String> netInterologFiles = new ArrayList<String>();
		for(String[] conf: interologConfig){
			System.err.println("Processing "+species_prefix+" "+conf[0]);
			netInterologFiles.add(buildInterologNetwork(workingDir, conf[0], workingDir+conf[1], workingDir+conf[2], workingDir+conf[3], workingDir+conf[4], species_prefix));
		}
		System.err.println("Merging files...");
		String allInterologFile = mergeInterolog(workingDir, netInterologFiles, species_prefix);
		networksByType.add(allInterologFile);
		networksByDb.addAll(netInterologFiles);
		System.err.println("Done.");
		
		//Complete network
		System.err.println("Creating complete network");
		String complete = buildCompleteNetwork(workingDir, allInterologFile, allDdiFile, species_prefix);
		
		final List<String> allNetworks = new ArrayList<String>();
		allNetworks.addAll(networksByDb);
		allNetworks.addAll(networksByType);
		allNetworks.add(complete);
		
		//Statistics for the network using gene ontology
		System.err.println("Running statistics");
		goCompartmentStats(workingDir, allNetworks, species_prefix);
		goInteractomeValidation(workingDir, allNetworks, species_prefix);
		System.err.println("Finished.");
		

	}
	
	/**
	 * Build interolog network for particular orthologous protein set
	 * 
	 * @throws Exception
	 */
	public static final String buildInterologNetwork(String workingDir, String speciesSource, String orthologFile, String swissprotMappingFile, String tremblMapplingFile, String intactFile, String speciesPrefix)throws Exception{
		SetMultimap<String, String> om = DelimitedIO.readEnsemblMap(orthologFile, 1, 0, null, 2);
		System.err.println("Orthologs: "+om.size());
		SetMultimap<String, String> um = DelimitedIO.readEnsemblMap(swissprotMappingFile, 1, 0, null, 2);
		SetMultimap<String, String> u2m = DelimitedIO.readEnsemblMap(tremblMapplingFile, 1, 0, null, 2);
		um.putAll(u2m);
		System.err.println("Uniprot mappings: "+um.size());
		
		
		BufferedReader br = new BufferedReader(new FileReader(intactFile));
		String line = br.readLine();
		Set<String> nodes = new HashSet<String>();
		Map<String,Double> edges = new HashMap<String, Double>();
		HashMultiset<String> times = HashMultiset.create();
		Pattern num = Pattern.compile("intact-miscore:(.*)");
		
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
					Matcher m = num.matcher(data[14]);
					if(m.find()){
						conf2 = Double.valueOf(m.group(1));
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
		
		String fileName = workingDir+speciesPrefix+"_interactome["+speciesSource+"].tab";
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		for(Entry<String, Double> ent: edges.entrySet()){
			bw.write(ent.getKey()+"\t"+ent.getValue().toString()+"\t"+String.valueOf(times.count(ent.getKey()))+"\n");
		}
		bw.flush();
		bw.close();
		return fileName;
	}
	
	public static final List<String> getAll(Collection<String> col, Multimap<String, String> map){
		Set<String> set = new HashSet<String>();
		for(String c : col){
			set.addAll(map.get(c));
		}
		return new ArrayList<String>(set);
	}

	/**
	 * Build domain-domain interaction network for individual database
	 * @throws Exception
	 */
	public static String buildDdiNetwork(String workingDir, String domainsFile, String interactingDomainsFile, String domainDbsuffix, String speciesPrefix, SetMultimap<String, ArrayKey<String>> map3) throws Exception{
		List<String[]> pairs = new ArrayList<String[]>();
		BufferedReader br = new BufferedReader(new FileReader(interactingDomainsFile));
		String line;// = br.readLine();
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			pairs.add(new String[]{data[0], data[1]});
			
		}
		br.close();
		
		br = new BufferedReader(new FileReader(domainsFile));
		SetMultimap<String, String> map = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			String dom = data[1].split("\\.")[0];
			map.put(dom, data[0].trim());
		}
		br.close();

		HashMultiset<String> ms = HashMultiset.create();
		SetMultimap<String, ArrayKey<String>> map2 = Multimaps.newSetMultimap(Maps.<String, Collection<ArrayKey<String>>> newHashMap(), DefaultSuppliers.<ArrayKey<String>> set());
		for(String [] pair : pairs){
			//System.err.println(pair[0]+"\t"+pair[1]);
			Collection<String> as = map.get(pair[0]);
			Collection<String> bs = map.get(pair[1]);
			for(String a :as){
				for(String b : bs){
					if(!a.equals(b)){
						if(a.compareTo(b) > 0){
							ms.add(a+"\t"+b);
							map2.put(a+"\t"+b, new ArrayKey<String>(new String[] {pair[0], pair[1]}));
							map3.put(a+"\t"+b, new ArrayKey<String>(new String[] {pair[0], pair[1]}));
						}
						else{
							ms.add(b+"\t"+a);
							map2.put(b+"\t"+a, new ArrayKey<String>(new String[] {pair[1], pair[0]}));
							map3.put(a+"\t"+b, new ArrayKey<String>(new String[] {pair[0], pair[1]}));
						}
					}
				}
			}
		}
		String fileName = workingDir+"\\"+speciesPrefix+"_ddi_["+domainDbsuffix+"].tab";
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		for(com.google.common.collect.Multiset.Entry<String> en : ms.entrySet()){
			bw.write(en.getElement()+"\t"+String.valueOf(en.getCount()+"\n"));
		}
		bw.flush();
		bw.close();
		
		String fileName2 = workingDir+"\\"+speciesPrefix+"_ddi_full_["+domainDbsuffix+"].tab";
		bw = new BufferedWriter(new FileWriter(fileName2));
		for(Entry<String, Collection<ArrayKey<String>>> en : map2.asMap().entrySet()){
			bw.write(en.getKey()+"\t"+String.valueOf(en.getValue().toString()+"\n"));
		}
		bw.flush();
		bw.close();
		return fileName2;
	}

	
	public static String buildCompleteNetwork(String workingDir, final String interologNetFile, final String ddinetFile, String species_prefix) throws Exception{
		Pattern p = Pattern.compile("\\[(.*?)\\]");
		List<String> files = new ArrayList<String>(){{add(interologNetFile);add(ddinetFile);}};
		List<String> ids = new ArrayList<String>();
		for(String f : files){
			Matcher m = p.matcher(f);
			m.find();
			ids.add(m.group(1));
		}
		Map<ArrayKey<String>, String[]> map = new HashMap<ArrayKey<String>, String[]>();
		BufferedReader br = new BufferedReader(new FileReader(files.get(0)));
		String line;
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			if(species_prefix.equals("Aspergillus_fumigatus")){
				data[0] = data[0].replace("CADAFUAT", "CADAFUAP");
				data[1] = data[1].replace("CADAFUAT", "CADAFUAP");
			}
			ArrayKey<String> key = ArrayKey.createSortedKey(data[0], data[1]);
			String [] value = new String [7];
			Arrays.fill(value, "");
			System.arraycopy(data, 0, value, 0, data.length);
			value[value.length-1] = ids.get(0);
			if(value[2].equals("NaN")){
				value[2]="";
			}
			map.put(key, value);
		}
		br.close();

		for(int i = 1; i < files.size(); i++){
			br = new BufferedReader(new FileReader(files.get(i)));
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				ArrayKey<String> key = ArrayKey.createSortedKey(data[0], data[1]);
				String [] d1 = map.get(key); 
				if(d1 == null){
					String [] value = new String [7];
					Arrays.fill(value, "");
					value[0] = data[0];
					value[1] = data[1];
					value[5] = data[2];
					value[value.length-1] = ids.get(i);
					map.put(key, value);
				}
				else{
					d1[d1.length-1] = "interologs,ddi";
					d1[5] = data[2];
					map.put(key, d1);
				}
			}
			br.close();
		}
		String fileName = workingDir+species_prefix+"_interactome[complete].tab";
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		bw.write("id1\tid2\tscore\tinteractions\tcount\tddi_dbs\tinference_method\n");
		for(String [] str : map.values()){
			bw.write(str[0]);
			for(int i = 1; i < str.length; i++){
				bw.write("\t"+str[i]);
			}
			bw.write("\n");
		}
		bw.flush();
		bw.close();
		return fileName;
	}
	
	public static String mergeDdi(String workingDir, List<String> files, String species_prefix) throws Exception{
		Pattern p = Pattern.compile("\\[(.*?)\\]");

		List<String> ids = new ArrayList<String>();
		for(String f : files){
			Matcher m = p.matcher(f);
			m.find();
			ids.add(m.group(1));
		}
		System.err.println(ids);

		Map<ArrayKey<String>, String[]> map = new HashMap<ArrayKey<String>, String[]>();
		BufferedReader br = new BufferedReader(new FileReader(files.get(0)));
		String line;
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			ArrayKey<String> key = ArrayKey.createSortedKey(data[0], data[1]);
			data[data.length-1] = ids.get(0);
			map.put(key, data);
			
		}
		br.close();
		
		for(int i = 1; i < files.size(); i++){
			br = new BufferedReader(new FileReader(files.get(i)));
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				ArrayKey<String> key = ArrayKey.createSortedKey(data[0], data[1]);
				String [] d1 = map.get(key); 
				if(d1 == null){
					data[data.length-1] = ids.get(i);
					map.put(key, data);
				}
				else{
					data[data.length-1] = ids.get(i)+","+d1[data.length-1];
					map.put(key, data);
				}
			}
			br.close();
		}
		
		String fileName = workingDir+species_prefix+"_interactome[ddi_all].tab";
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		for(String [] str : map.values()){
			bw.write(str[0]);
			for(int i = 1; i < str.length; i++){
				bw.write("\t"+str[i]);
			}
			bw.write("\n");
		}
		bw.flush();
		bw.close();
		return fileName;
	}

	public static String mergeInterolog(String workingDir, List<String> files, String species_prefix) throws Exception{
		Pattern p = Pattern.compile("\\[(.*?)\\]");

		List<String> ids = new ArrayList<String>();
		for(String f : files){
			Matcher m = p.matcher(f);
			m.find();
			ids.add(m.group(1));
		}
		Map<ArrayKey<String>, String[]> map = new HashMap<ArrayKey<String>, String[]>();
		BufferedReader br = new BufferedReader(new FileReader(files.get(0)));
		String line;
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			ArrayKey<String> key = ArrayKey.createSortedKey(data[0], data[1]);
			data = Arrays.copyOf(data, data.length+1);
			data[data.length-1] = ids.get(0);
			map.put(key, data);
		}
		br.close();
		
		for(int i = 1; i < files.size(); i++){
			br = new BufferedReader(new FileReader(files.get(i)));
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				ArrayKey<String> key = ArrayKey.createSortedKey(data[0], data[1]);
				String [] d1 = map.get(key); 
				if(d1 == null){
					data = Arrays.copyOf(data, data.length+1);
					data[data.length-1] = ids.get(i);
					map.put(key, data);
				}
				else{
					data = Arrays.copyOf(data, data.length+1);
					data[data.length-1] = ids.get(i);
					data[data.length-1] = data[data.length-1]+","+d1[data.length-1];
					Double dScore = Double.valueOf(data[2]);
					Double d1Score =  Double.valueOf(d1[2]);
					if(Double.isNaN(dScore) && !Double.isNaN(d1Score)){
						data[2] = d1Score.toString();
					}
					else if(!Double.isNaN(dScore) && Double.isNaN(d1Score)){
						data[2] = dScore.toString();
					}
					else if(dScore > d1Score){
						data[2] = dScore.toString();
					}
					else{
						data[2] = d1Score.toString();
					}
					try{
						data[3] = String.valueOf(Double.valueOf(data[3]) + Double.valueOf(d1[3]));
					}
					catch(NumberFormatException e){
						System.err.println(Arrays.asList(data).toString() +"\t"+data[3]+"\t"+d1[3]);
					}
					map.put(key, data);
					
					
				}
			}
			br.close();
		}
		String fileName = workingDir+species_prefix+"_interactome[interalogs_all].tab";
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		for(String [] str : map.values()){
			bw.write(str[0]);
			for(int i = 1; i < str.length; i++){
				bw.write("\t"+str[i]);
			}
			bw.write("\n");
		}
		bw.flush();
		bw.close();
		return fileName;
	}

	
	public static void goCompartmentStats(String workingDir, List<String> networks, String species_prefix) throws Exception{
		Map<String, String> map = new HashMap<String, String>(){{
			put("extracellular region", "GO:0005576");
			put("cytoplasm", "GO:0005737");
			put("nucleus", "GO:0005634");
			put("mitochondrion", "GO:0005739");
			put("endoplasmic reticulum", "GO:0005783");
			put("Golgi apparatus", "GO:0005794");
			put("fungal-type vacuole", "GO:0000324");
			put("fungal-type cell wall", "GO:0009277");
		}};
		
		Pattern p = Pattern.compile("\\[(.*?)\\]");

		
		List<String> net_names = new ArrayList<String>();
		for(String f : networks){
			Matcher m = p.matcher(f);
			m.find();
			net_names.add(m.group(1));
		}

		OntologyGraph onto = new OntologyGraph();
		String obo = workingDir+"go\\go.obo";
		String [] rels = new String[]{"is_a"};
		String anno = workingDir+"go\\goa.tab";
		onto.readObo(obo, "cellular_component"); //"biological_process", "molecular_function", "cellular_component"
		onto.indexHierarchy(rels);
		onto.loadAnnotationsFromTab(anno);
		StringBuffer sb = new StringBuffer();
		
		SetMultimap<String, Term> ccs = Multimaps.newSetMultimap(Maps.<String, Collection<Term>> newHashMap(), DefaultSuppliers.<Term> set());
		
		for(Entry<String, String> ent : map.entrySet()){
			Term t =  onto.getTerm(ent.getValue());
			int ans = onto.getAnnotatedEntities(t).size();
			for(Entity e : onto.getAnnotatedEntities(t)){
				ccs.put(e.id.toString(), t);
			}
			sb.append(t.name+"\t"+String.valueOf(ans)+"\n");
		}
		List<String> ents = new ArrayList<String>(ccs.keySet());
		sb.append(ents.size()+"\n");
	
		double total = 0;
		double same = 0;
		for(int i = 0; i < networks.size(); i++){
			BufferedReader br = new BufferedReader(new FileReader(networks.get(i)));
			String line;
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				Collection<Term> as = ccs.get(data[0]);
				Collection<Term> bs = ccs.get(data[1]);
				if(!as.isEmpty() && !bs.isEmpty()){
					total+= 1d;
					if(hasMatch(as, bs)){
						same+=1d;
					}
				}
			}
			br.close();
			sb.append(net_names.get(i)+"\t"+String.valueOf(same/total)+"\t"+String.valueOf(same)+"\t"+String.valueOf(total)+"\n");
		}
		DelimitedIO.saveToFile(sb, workingDir+species_prefix+"_compartment_stats.txt");
	}
		
	private static class MicaEvaluator{
		private Map<String, OntologyGraph> map = new HashMap<String, OntologyGraph>();
		Map<String, Map<String, Set<Term>>> cache = new HashMap<String, Map<String, Set<Term>>>();
		
		public MicaEvaluator(String workingDir, String species_prefix) throws Exception{
			String [] aspects = new String[]{"biological_process", "molecular_function"}; //, "cellular_component"
			String obo = workingDir+"go\\go.obo";
			String anno = workingDir+"goa\\"+species_prefix+"_goa.tab";
			for(String aspect : aspects){
				cache.put(aspect, new HashMap<String, Set<Term>>());
				OntologyGraph onto = new OntologyGraph();
				String [] rels = new String[]{"is_a", "part_of"};
				onto.readObo(obo, aspect); //"biological_process", "molecular_function", "cellular_component"
				onto.indexHierarchy(rels);
				onto.loadAnnotationsFromTab(anno);
				onto.calculateIc();
				
				map.put(aspect, onto);
			}
		}

		public Double getMicaScore(String aspect, String a, String b){
			OntologyGraph onto = map.get(aspect);
			Entity e1 = onto.convert(a);
			Entity e2 = onto.convert(b);
			
			if(e1 == null || e2 == null){
				return null;
			}
			
			Map<String, Set<Term>> tmp = cache.get(aspect);
			
			Set<Term> t1 = tmp.get(a);
			if(t1 == null){
				t1 = onto.getAllAnnotationTerms(e1);
				tmp.put(a, t1);
			}
			
			Set<Term> t2 = tmp.get(b);
			if(t2 == null){
				t2 = onto.getAllAnnotationTerms(e2);
				tmp.put(b, t2);
			}
			Term m = onto.getMaxMica(t1, t2);
			if(m == null){
				return 0d;
			}
			return m.getIc();
		}
		
		public Set<String> getAllAnnotated(String aspect){
			OntologyGraph onto = map.get(aspect);
			Set<String> set = new HashSet<String>(onto.entityToId(onto.getAnnotations().keySet()));
			return set;
		}
	}

	private static class CompartmentEvaluator{
		private SetMultimap<String, Term> ccs = Multimaps.newSetMultimap(Maps.<String, Collection<Term>> newHashMap(), DefaultSuppliers.<Term> set());
		private OntologyGraph onto = new OntologyGraph();
		
		public CompartmentEvaluator(String workingDir, String species_prefix) throws Exception{
			Map<String, String> map = new HashMap<String, String>(){{
				put("extracellular region", "GO:0005576");
				put("cytoplasm", "GO:0005737");
				put("nucleus", "GO:0005634");
				put("mitochondrion", "GO:0005739");
				put("endoplasmic reticulum", "GO:0005783");
				put("Golgi apparatus", "GO:0005794");
				put("fungal-type vacuole", "GO:0000324");
				put("fungal-type cell wall", "GO:0009277");
			}};

			OntologyGraph onto = new OntologyGraph();
			String obo = workingDir+"go\\go.obo";
			String [] rels = new String[]{"is_a"};
			String anno = workingDir+"goa\\"+species_prefix+"_goa.tab";
			onto.readObo(obo, "cellular_component"); //"biological_process", "molecular_function", "cellular_component"
			onto.indexHierarchy(rels);
			onto.loadAnnotationsFromTab(anno);

			for(Entry<String, String> ent : map.entrySet()){
				Term t =  onto.getTerm(ent.getValue());
				for(Entity e : onto.getAnnotatedEntities(t)){
					ccs.put(e.id.toString(), t);
				}
			}
		}
		
		public Integer eval(String a, String b){
			Collection<Term> as = ccs.get(a);
			Collection<Term> bs = ccs.get(b);
			if(!as.isEmpty() && !bs.isEmpty()){
				if(hasMatch(as, bs)){
					return 1;
				}
				return 0;
			}
			return null;
		}
		
		public Set<String> getAnnotated(){
			Set<String> set = new HashSet<String>(ccs.keySet());
			return set;
		}
	}
	
	public static boolean hasMatch(Collection as, Collection bs){
		for(Object a : as){
			if(bs.contains(a)){
				return true;
			}
		}
		return false;
	}
	
	public static void goInteractomeValidation(String workingDir, List<String> networks, String species_prefix) throws Exception{
		Pattern p = Pattern.compile("\\[(.*?)\\]");
		List<String> net_names = new ArrayList<String>();
		for(String f : networks){
			Matcher m = p.matcher(f);
			m.find();
			net_names.add(m.group(1));
		}
		System.err.println(net_names);
		StringBuffer sb = new StringBuffer();
		String [] aspects = new String[]{"biological_process", "molecular_function", "cellular_component"};
		String obo = workingDir+"go\\go.obo";
		String anno = workingDir+"go\\goa.tab";
		for(String aspect : aspects){
			OntologyGraph onto = new OntologyGraph();
			String [] rels = new String[]{"is_a", "part_of"};
			onto.readObo(obo, aspect); //"biological_process", "molecular_function", "cellular_component"
			onto.indexHierarchy(rels);
			onto.loadAnnotationsFromTab(anno);
			onto.calculateIc();
			for(int z = 0; z < net_names.size(); z++){
				Map<ArrayKey<Entity>, Double> annotatedEdges = new HashMap<ArrayKey<Entity>, Double>();
				String net = networks.get(z);
				
				String netName = net_names.get(z);
				System.err.println("Running "+species_prefix+" "+netName+" network on "+aspect+" -- file: "+networks.get(z));
				BufferedReader br = new BufferedReader(new FileReader(net));
				String line;
				List<Double> real = new ArrayList<Double>();
				
				BufferedWriter bw = new BufferedWriter(new FileWriter(workingDir+species_prefix+"_interactome_validataion_"+aspect+".tab", true));

				List<Entity> ents = new ArrayList<Entity>();
				int skip1 = 0;
				int skip2 = 0;

				while ((line = br.readLine()) != null) {
					String[] data = line.split("\t");
					Entity e1 = onto.convert(data[0]);
					Entity e2 = onto.convert(data[1]);
					//if("D:\\projects\\path_net_combined\\fg_interactome[ddi_all].tab".equals(net)){
						//System.err.println(data[0]+"="+(e1 == null)+" -- "+data[1]+"="+(e1 == null)+" >> "+line);
					//}
					if(e1 == null || e2 == null){
						continue;
					}
					skip1++;
					ArrayKey<Entity> key = ArrayKey.createSortedKey(e1, e2);
					Double d = annotatedEdges.get(key);

					if(d == null){
						Set<Term> t1 = onto.getAllAnnotationTerms(e1);
						if( !t1.isEmpty()){
							ents.add(e1);
						}

						Set<Term> t2 = onto.getAllAnnotationTerms(e2);
						if(!t2.isEmpty()){
							ents.add(e2);
						}
						
						if(t1.isEmpty() || t2.isEmpty()){
							continue;
						}
						skip2++;
						Term m = onto.getMaxMica(t1, t2);
						if(m != null){
							real.add(m.getIc());
							bw.write(data[0]+"_"+data[1]+"\t"+netName+"\t"+aspect+"\t"+String.valueOf(m.getIc())+"\n");
							annotatedEdges.put(key, m.getIc());
						}
						else{
							real.add(0.0d);
							bw.write(data[0]+"_"+data[1]+"\t"+netName+"\t"+aspect+"\t"+String.valueOf(0.0d)+"\n");
							annotatedEdges.put(key, 0.0d);
						}
					}
					else{
						bw.write(data[0]+"_"+data[1]+"\t"+netName+"\t"+aspect+"\t"+String.valueOf(d)+"\n");
					}

				}
				br.close();
				bw.flush();
				bw.close();

				List<Double> random = new ArrayList<Double>();
				Random r = new Random(System.currentTimeMillis());
				while(random.size() < real.size()){
					int i = r.nextInt(ents.size());
					int j = r.nextInt(ents.size());
					if(i == j){
						continue;
					}
					ArrayKey<Entity> key = ArrayKey.createSortedKey(ents.get(i), ents.get(j));
					Double d = annotatedEdges.get(key);
					if(d == null){
						Set<Term> t1 = onto.getAllAnnotationTerms(ents.get(i));
						Set<Term> t2 = onto.getAllAnnotationTerms(ents.get(j));
						Term m = onto.getMaxMica(t1, t2);
						if(m != null){
							random.add(m.getIc());
							annotatedEdges.put(key, m.getIc());
						}
						else{
							random.add(0.0d);
							annotatedEdges.put(key, 0.0d);
						}
					}
				}
				System.err.println("Random size: "+random.size());
				double [] dreal = Doubles.toArray(real);
				double [] drandom = Doubles.toArray(random);
				sb.append(netName+" - "+aspect+"\n");
				sb.append("real:\t"+StatUtils.mean(dreal)+"\n");
				sb.append("random:\t"+StatUtils.mean(drandom)+"\n");
				MannWhitneyUTest mwt = new MannWhitneyUTest();
				sb.append("p-value:\t"+mwt.mannWhitneyUTest(dreal, drandom)+"\t"+mwt.mannWhitneyU(dreal, drandom)+"\n");
				sb.append("\n");
			}
		}
		DelimitedIO.saveToFile(sb, workingDir+species_prefix+"_interactome_validataion_summary.txt");
	}
}
