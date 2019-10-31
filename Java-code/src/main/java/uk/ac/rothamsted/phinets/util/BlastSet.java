package uk.ac.rothamsted.phinets.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Multimaps;
import com.google.common.collect.SetMultimap;

public class BlastSet {
	
	private final Map<String, String> id2subset = new HashMap<String, String>();
	private final SetMultimap<String, String> subset2id = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
	private final Map<ArrayKey<String>,Hit> hits = new HashMap<ArrayKey<String>,Hit>();
	private final Map<String, Map<ArrayKey<String>,Hit>> index = new HashMap<String, Map<ArrayKey<String>,Hit>>();
	
	
	
	public void generateScores(String outDir, double min) throws Exception{
		Map<String, BufferedWriter> bws = new HashMap<String, BufferedWriter>();
		for(String s : subset2id.keySet()){
			BufferedWriter bw = new BufferedWriter(new FileWriter(outDir+s.split("\\.")[0]+"_seq.tab"));
			bws.put(s, bw);
			bw.write("species\tid\tinscore\toutscore\n");
		}
		ArrayKey<String> masterKey = new ArrayKey<String>();
		SetMultimap<String, String> processed = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		for(Entry<String, Map<ArrayKey<String>,Hit>> ent1 : index.entrySet()){
			String id = ent1.getKey();
			String subset = id2subset.get(id);
			BufferedWriter bw = bws.get(subset);
			double bestIn = Double.MIN_VALUE;
			Hit hIn = null;
			double bestOut = Double.MIN_VALUE;
			Hit hOut = null;
			
			for(Entry<ArrayKey<String>,Hit> ent : ent1.getValue().entrySet()){
				String [] arr = ent.getKey().getArray();
				if(!arr[0].equals(id)){
					continue;
				}
				Hit a = ent.getValue();
				double ea = a.evalue==0?min:a.evalue;
				String [] temp = new String[]{arr[1],arr[0]};
				masterKey.setArray(temp);
				if(a.subsetTo.equals(subset) && a.subsetFrom.equals(subset)){
					Hit b = ent1.getValue().get(masterKey);
					if(b!= null){
						double eb = b.evalue==0?min:b.evalue;
						double score = (-Math.log10(ea)-Math.log10(eb))/2d;
						if(bestIn < score){
							bestIn = score;
							//hIn = a;
						}
					}
				}
				else{ 
					Hit b = ent1.getValue().get(masterKey);
					if(b!= null){
						double eb = b.evalue==0?min:b.evalue;
						double score = (-Math.log10(ea)-Math.log10(eb))/2d;
						if(bestOut < score){
							bestOut = score;
							//hOut = a;
						}
					}
				}
			}
			processed.put(subset, id);
			bw.write(subset);
			bw.write("\t");
			bw.write(id);
			bw.write("\t");
			if(bestIn == Double.MIN_VALUE){
				bw.write("5.0");
			}
			else{
				bw.write(String.valueOf(bestIn));
				//bw.write("\t"+hIn.from+"--"+hIn.to);
			}
			bw.write("\t");
			if(bestOut == Double.MIN_VALUE){
				bw.write("5.0");
			}
			else{
				bw.write(String.valueOf(bestOut));
				//bw.write("\t"+hOut.from+"--"+hOut.to);
			}
			bw.write("\n");
			bw.flush();
		}
		SetMultimap<String, String> all = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		all.putAll(subset2id);
		for(Entry<String, BufferedWriter> ent: bws.entrySet()){
			BufferedWriter bw = ent.getValue();
			Set<String> set = all.get(ent.getKey());
			set.removeAll(processed.get(ent.getKey()));
			for(String s : set){
				bw.write(ent.getKey());
				bw.write("\t");
				bw.write(s);
				bw.write("\t");
				bw.write("5.0");
				bw.write("\t");
				bw.write("5.0");
				bw.write("\n");
			}
			bw.flush();
			bw.close();
		}

	}
	
	public void buildSubsetIndex(String faDir) throws Exception{
		for(File f : new File(faDir).listFiles()){
			if(f.isDirectory()){
				continue;
			}
			String name = f.getName();
			BufferedReader br = new BufferedReader(new FileReader(f));
			String line;
			while ((line = br.readLine()) != null) {
				if(!line.startsWith(">")){
					continue;
				}
				String[] data = line.split(" ");
				data[0]=data[0].replace(">", "");
				//System.err.println(data[0]);
				id2subset.put(data[0], name);
				subset2id.put(name, data[0]);
			}
			br.close();
			System.err.println(name+"\t"+subset2id.get(name).size());
		}
	}
	
	public double loadBlast(String file) throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		double min = Double.MAX_VALUE;
		while ((line = br.readLine()) != null) {
			new Hit(line, id2subset, hits);
			String [] data = line.split("\t");
			double temp = Double.valueOf(data[8]);
			if(temp > 0 && temp < min){
				min = temp;
			}
			
		}
		br.close();
		System.err.println(min);
		System.err.println(hits.size());
		
		for(Entry<ArrayKey<String>,Hit> ent :hits.entrySet()){
			Map<ArrayKey<String>,Hit> m2 = index.get(ent.getKey().getArray()[0]);
			if(m2 == null){
				m2 = new HashMap<ArrayKey<String>,Hit>();
				index.put(ent.getKey().getArray()[0], m2);
			}
			m2.put(ent.getKey(), ent.getValue());
			Map<ArrayKey<String>,Hit> m1 = index.get(ent.getKey().getArray()[1]);
			if(m1 == null){
				m1 = new HashMap<ArrayKey<String>,Hit>();
				index.put(ent.getKey().getArray()[1], m1);
			}
			m1.put(ent.getKey(), ent.getValue());
		}
		System.err.println(index.size());
		System.err.println("Done.");
		return min;
	}
	
	private static class Hit {
		public final String from;
		public final String to;
		public final double evalue;
		public final double bitscore;
		public final double fromLength;
		public final double toLength;
		public final double alignmentLength;
		public final String subsetFrom;
		public final String subsetTo;
		
		public Hit(String line, Map<String, String> id2subset, Map<ArrayKey<String>, Hit> hits){
			String [] data = line.split("\t");
			from = data[1].split("\\|")[0];
			to = data[4].split("\\|")[0];
			ArrayKey<String> key = new ArrayKey<String>(new String[]{from, to});
			evalue = Double.valueOf(data[8]);
			bitscore = Double.valueOf(data[7]);
			fromLength = Double.valueOf(data[2]);
			toLength = Double.valueOf(data[5]);
			alignmentLength = Double.valueOf(data[9]);
			subsetFrom = id2subset.get(from);
			//if(subsetFrom == null){
			//	System.err.println(from+" :: "+to+" :: "+line);
			//}
			subsetTo = id2subset.get(to);
			
			double cov = alignmentLength/(Math.min(fromLength, toLength));

			if(cov > 0.3 && !from.equals(to)){
				Hit h = hits.get(key);
				if(h == null || evalue < h.evalue){
					hits.put(key, this);
				}
			}
		}
	}
}
