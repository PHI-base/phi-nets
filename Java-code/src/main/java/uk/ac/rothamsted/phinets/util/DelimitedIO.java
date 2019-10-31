package uk.ac.rothamsted.phinets.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.collect.SetMultimap;

import uk.ac.rothamsted.phinets.util.ArrayKey;

public class DelimitedIO {

	public static void removeBlanks(String fin, String fout, boolean removeHeader)throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(fin));
		if(removeHeader){
			br.readLine();
		}
		String line;
		BufferedWriter bw = new BufferedWriter(new FileWriter(fout));
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
			if(data.length < 2 || data[0].trim().isEmpty() || data[1].trim().isEmpty()){
				continue;
			}
			bw.write(line+"\n");
		}
		br.close();
		bw.flush();
		bw.close();
	}
	
	public static void saveMapDiamondFormat(SetMultimap<String, String> map, String file, int minValues) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		bw.write("#Key values\n");
		for(Entry<String, Collection<String>> ent : map.asMap().entrySet()){
			if(ent.getValue().size() < minValues){
				continue;
			}
			bw.write(ent.getKey());
			bw.write("  ");
			boolean first = true;
			for(String s : ent.getValue()){
				if(!first){
					bw.write("/");
				} 
				first = false;
				bw.write(s);
			}
			bw.write("\n");
		}
		bw.flush();
		bw.close();
	}

	//List<Entry<String, Double>>

	public static void saveToFile(List<Entry<String, Double>> li, String file, String delimiter) throws Exception {
		DataOutputStream out = new DataOutputStream(new FileOutputStream(file));
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out));
		for (Entry ent : li) {
			bw.write(ent.getKey().toString());
			bw.write("\t");
			bw.write(ent.getValue().toString());
			bw.write("\n");
		}
		bw.flush();
		out.close();
	}

	public static SetMultimap<String, String> readNetwork(String file, String  delimiter)throws Exception{
		SetMultimap<String, String> map = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		int edges = 0;
		while((strLine = br.readLine()) != null){
			String [] data = strLine.split(delimiter);
			if(data.length > 1 && data.length  > 0){
				String a = data[0].trim();
				String b = data[1].trim();
				if(!a.equals("") && !b.equals("") && !a.equals(b)){
					map.put(a, b);
					map.put(b, a);
					edges++;
				}
			}
		}
		System.err.println("Network with "+edges+" edges");
		in.close();
		return map;
	}
	
	public static SetMultimap<String, String> readNetwork(String file, String  delimiter, boolean directed)throws Exception{
		SetMultimap<String, String> map = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		int edges = 0;
		while((strLine = br.readLine()) != null){
			String [] data = strLine.split(delimiter);
			if(data.length > 1 && data.length  > 0){
				String a = data[0].trim();
				String b = data[1].trim();
				if(!a.equals("") && !b.equals("")){
					map.put(a, b);
					if(!directed){
						map.put(b, a);
					}
					edges++;
				}
			}
		}
		System.err.println("Network with "+edges+" edges");
		in.close();
		return map;
	}
	
	public static Set<ArrayKey<String>> toArrayKeyList(SetMultimap<String, String> ... maps){
		Set<ArrayKey<String>> list = new HashSet<ArrayKey<String>>();
		for(SetMultimap<String, String> map : maps){
			for(Entry<String, Collection<String>> ent : map.asMap().entrySet()){
				for(String s : ent.getValue()){
					list.add(new ArrayKey<String>(true, ent.getKey(), s));
				}
			}
		}
		return list;
	}
	
	public static void save(Collection<ArrayKey<String>> list, String file) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		for(ArrayKey key : list){
			Object [] o = key.getArray();
			if(o.length == 0){
				continue;
			}
			bw.write(o[0].toString());
			for(int i =1; i< o.length; i++){
				bw.write("\t");
				bw.write(o[i].toString());
			}
			bw.write("\n");
		}
		bw.flush();
		bw.close();

	}
	
	public static <V> ListMultimap<Integer, V> sortIntoBinsInt(ListMultimap<Integer, V> input, int numBins, Integer maxLimit, boolean discard){
		ListMultimap<Integer, V> list = ArrayListMultimap.create();
		int max = maxLimit==null?Collections.max(input.keySet()):maxLimit;
		int min = Collections.max(input.keySet());
		final int binSize = (max - min)/numBins;
		List<Integer> current = new ArrayList<Integer>(input.keySet());
		Collections.sort(current);
		int start = min;
		
		for(Integer i : current){
			while(i != null){
				if(i <= start+binSize){
					list.putAll(start+(binSize/2), input.get(i));
					i = null;
				}
				else if(i > max){
					if(!discard){
						list.putAll(start+(binSize/2), input.get(i));
					}
					i = null;
				}
				else{
					start += binSize;
				}
			}
		}
		return list;
	}
	
	public static <V> ListMultimap<Double, V> sortIntoBinsDouble(ListMultimap<Double, V> input, int numBins, Double maxLimit, boolean discard){
		ListMultimap<Double, V> list = ArrayListMultimap.create();
		double max = maxLimit==null?Collections.max(input.keySet()):maxLimit;
		double min = Collections.max(input.keySet());
		final double binSize = (max - min)/((double)numBins);
		List<Double> current = new ArrayList<Double>(input.keySet());
		Collections.sort(current);
		double start = min;
		
		for(Double i : current){
			while(i != null){
				if(i <= start+binSize){
					list.putAll(start+(binSize/2d), input.get(i));
					i = null;
				}
				else if(i > max){
					if(!discard){
						list.putAll(start+(binSize/2d), input.get(i));
					}
					i = null;
				}
				else{
					start += binSize;
				}
			}
		}
		return list;
	}
	
	public static void mergeColumns(String dir, String out, String empty) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		List<BufferedReader> brs = new ArrayList<BufferedReader>();
		List<String> headers = new ArrayList<String>();
		for(File f : new File(dir).listFiles()){
			if(f.isDirectory()){
				continue;
			}
			BufferedReader br = new BufferedReader(new FileReader(f));
			brs.add(br);
			headers.add(f.getName().replaceAll("\\.", "_"));
		}
		
		for(int i = 0; i < brs.size(); i++){
			if(i > 0){
				bw.write("\t");
			}
			bw.write(headers.get(i));
		}
		bw.write("\n");
		
		String [] data = new String[brs.size()];
		boolean done = false;
		while(!done){
			for(int i = 0; i < brs.size(); i++){
				done = true;
				data[i] = brs.get(i).readLine();
				if(data[i] != null){
					done = false;
				}
			}
			if(!done){
				for(int i = 0; i < brs.size(); i++){
					if(i > 0){
						bw.write("\t");
					}
					if(data[i] != null){
						bw.write(data[i]);
					}
					else if(empty != null){
						bw.write(empty);
					}
				}
				bw.write("\n");
			}
		}
		bw.flush();
		bw.close();

		for(BufferedReader br : brs){
			br.close();
		}
	}
	
	public static Set<String> toSet(String in, String split){
		Set<String> set = new HashSet<String>();
		String [] data = in.split(split);
		for(String d : data){
			d = d.trim();
			if(!d.isEmpty()){
				set.add(d);
			}
		}
		return set;
	}
	
	public static String makeString(Collection<String> set, String sep){
		 StringBuffer sb = new StringBuffer();
		 boolean first = true;
		 for(String s : set){
			 if(!first){
				 sb.append(sep);
			 }
			 first = false;
			 sb.append(s);
		 }
		 return sb.toString();
	}
	
	public static void padMap(ListMultimap map, Object pad, Collection<Object> keys){
		Set<Object> set = new HashSet<Object>(keys);
		set.removeAll(map.keys());
		for(Object o: set){
			map.put(o, pad);
		}
		padMap(map, pad);
	}
	
	public static void padMap(ListMultimap<Object, Object> map, Object pad){
		int max = 0;
		for(Collection c : map.asMap().values()){
			if(max < c.size()){
				max = c.size();
			}
		}
		Set<Object> keys = new HashSet<Object>(map.keys());
		for(Object key : keys){
			for(int i = 0; i < (max - map.get(key).size()); i++){
				map.put(key, pad);
			}
		}
		
	}
	
	
	public static List<String> readColumn(String file, int column, String delimiter) throws Exception {
		List<String> list = new ArrayList<String>();
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			if(delimiter != null){
				String[] data = strLine.split(delimiter);
				if (data.length > column) {
					list.add(data[column]);
				}	
			}
			else{
				list.add(strLine);
			}
		}
		in.close();
		return list;
	}
	
	public static List<Integer> readIntColumn(String file, int column, String delimiter, boolean header) throws Exception {
		List<Integer> list = new ArrayList<Integer>();
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		if(header){
			br.readLine();
		}
		while ((strLine = br.readLine()) != null) {
			if(delimiter != null){
				String[] data = strLine.split(delimiter);
				if (data.length > column && !data[column].equals("")) {
					list.add(Integer.valueOf(data[column]));
				}	
			}
			else if(!strLine.equals("")){
				list.add(Integer.valueOf(strLine));
			}
		}
		in.close();
		return list;
	}
	
	public static Multimap<String, String> readSimpleMap(Multimap<String, String> map, String file, String  delimiter, int col1, int col2)throws Exception{
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while((strLine = br.readLine()) != null){
			String [] data = strLine.split(delimiter);
			if(data.length > col2 && data.length  > col1){
				String a = data[col1].trim();
				String b = data[col2].trim();
				if(!a.equals("") && !b.equals("")){
					map.put(a, b);
				}
			}
		}
		in.close();
		return map;
	}
	
	public static SetMultimap<String, String> readEnsemblMap(String file, int col1, int col2, String valueIfKeep, int filterCol)throws Exception{
		String  delimiter = "\t";
		SetMultimap<String, String> map = Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String> set());
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine =br.readLine();
		while((strLine = br.readLine()) != null){
			String [] data = strLine.split(delimiter);
			if(valueIfKeep != null && (data.length < filterCol || !data[filterCol].equals(valueIfKeep))){
				continue;
			}
			if(data.length > col2 && data.length  > col1){
				String a = data[col1].trim();
				String b = data[col2].trim();
				if(!a.equals("") && !b.equals("")){
					map.put(a, b);
				}
			}
		}
		in.close();
		return map;
	}
	
	//
	
	public static Map<ArrayKey<String>, String> readCompoundKeyMap(String file, String  delimiter, int value, boolean sort, int ... keys)throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(file));
		Map<ArrayKey<String>, String> map = new HashMap<ArrayKey<String>, String>();
		String strLine;
		int max = max(keys);
		while((strLine = br.readLine()) != null){
			String [] data = strLine.split(delimiter);
			if(data.length > value && data.length  > max){
				String v = data[value].trim();
				if(!v.equals("")){
					String [] key = new String[keys.length];
					for(int i = 0; i< keys.length ; i++){
						key[i] = data[keys[i]];
					}
					//System.err.println(Arrays.asList(key).toString());
					ArrayKey<String> ak = new ArrayKey<String>(key);
					if(sort){
						ak.sort();
					}
					map.put(ak, v);
				}
			}
		}
		br.close();
		return map;
	}
	
	public static Map<String, Integer[]> readCoordinatesMap(String file, String  delimiter, int key, int ... values)throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(file));
		Map<String, Integer[]> map = new HashMap<String, Integer[]>();
		String strLine;
		int max = max(values);
		while((strLine = br.readLine()) != null){
			String [] data = strLine.split(delimiter);
			if(data.length > key && data.length  > max && !data[key].equals("")){
				Integer [] val = new Integer[values.length];
				for(int i = 0; i < values.length; i++){
					try{
						val[i] = Integer.valueOf(data[values[i]]);
					}
					catch(NumberFormatException e){
						System.err.println("Skipped "+data[key]+" due to "+data[values[i]]);
						continue;
					}
				}
				map.put(data[key], val);
			}
		}
		br.close();
		return map;
	}
	
	public static Map<String, Map<String, Integer[]>> readSplitCoordinatesMap(String file, String  delimiter, int split, int key, int ... values)throws Exception{
		BufferedReader br = new BufferedReader(new FileReader(file));
		Map<String, Map<String, Integer[]>> result = new HashMap<String, Map<String, Integer[]>>();
		String strLine;
		int max = max(values);
		while((strLine = br.readLine()) != null){
			String [] data = strLine.split(delimiter);
			if(data.length > key && data.length  > max && !data[key].equals("")){
				Map<String, Integer[]> map = result.get(data[split]);
				if(map == null){
					map = new HashMap<String, Integer[]>();
					result.put(data[split], map);
				}
				Integer [] val = new Integer[values.length];
				for(int i = 0; i < values.length; i++){
					try{
						val[i] = Integer.valueOf(data[values[i]]);
					}
					catch(NumberFormatException e){
						System.err.println("Skipped "+data[key]+" due to "+data[values[i]]);
						continue;
					}
				}
				map.put(data[key], val);
			}
		}
		br.close();
		return result;
	}
	
	public static <V,S> Set<S> getValues(Collection<V> ids, Map<V, S> map){
		Set<S> set = new HashSet<S>();
		for(V id : ids){
			set.add(map.get(id));
		}
		return set;
	}
	
	public static <V> String makeList(Collection<V> col, String delimiter){
		StringBuffer sb = new StringBuffer();
		boolean first = true;
		for(V v : col){
			if(!first){
				sb.append(delimiter);
			}
			first = false;
			sb.append(v.toString());
		}
		return sb.toString();
	}
	
	public static<A,B> SetMultimap<B, A> invert(SetMultimap<A, B> map){
		SetMultimap<B, A> out = Multimaps.newSetMultimap(Maps.<B, Collection<A>> newHashMap(), DefaultSuppliers.<A> set());
		for (Entry<A, Collection<B>> ent : map.asMap().entrySet()) {
			for(B s : ent.getValue()){
				out.put(s, ent.getKey());
			}
		}
		return out;
		
	}

	public static int max(int [] arr){
		int max = Integer.MIN_VALUE;
		for(int a : arr){
			if(a > max){
				max = a;
			}
		}
		return max;
	}
	
	public static Map<String, String> readSimpleMap(Map<String, String> map, String file, String  delimiter, int col1, int col2)throws Exception{
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while((strLine = br.readLine()) != null){
			String [] data = strLine.split(delimiter);
			if(data.length > col2 && data.length  > col1){
				String a = data[col1].trim();
				String b = data[col2].trim();
				if(!a.equals("") && !b.equals("")){
					map.put(a, b);
				}
			}
		}
		in.close();
		return map;
	}
	
	public static StringBuffer summary(SetMultimap<String, String> map){
		StringBuffer sb = new StringBuffer();
		for(Entry<String,Collection<String>> ent : map.asMap().entrySet()){
			sb.append(ent.getKey());
			sb.append("\t");
			sb.append(ent.getValue().size());
			sb.append("\n");
		}
		return sb;
	}
	
	public static List<String> readColumn(InputStream file, int column, String delimiter) throws Exception {
		List<String> list = new ArrayList<String>();
		DataInputStream in = new DataInputStream(file);
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			if(delimiter != null){
				String[] data = strLine.split(delimiter);
				if (data.length > column) {
					list.add(data[column]);
				}	
			}
			else{
				list.add(strLine);
			}
		}
		in.close();
		return list;
	}

	public static void saveToFile(Collection<?> data, String file, String delimiter) throws Exception {
		DataOutputStream out = new DataOutputStream(new FileOutputStream(file));
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out));
		for (Object object : data) {
			if ((object instanceof Collection)) {
				writeRow(bw, (Collection) object, delimiter);
			} else {
				bw.write(object.toString());
				bw.flush();
			}
			bw.write("\n");
			bw.flush();
		}
		out.close();
	}

	private static void writeRow(BufferedWriter bw, Collection<?> row, String delimiter) throws Exception {
		if (row.size() == 0) {
			return;
		}
		Iterator<?> it = row.iterator();
		bw.write(it.next().toString());
		while (it.hasNext()) {
			bw.write(delimiter);
			bw.write(it.next().toString());
		}
		bw.flush();
	}
	
	public static <A, B> void saveVectors(ListMultimap<A, B> lm, String file, String prefix, String na) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		int max = 0;
		boolean first = true;
		for (Entry<A, Collection<B>> ent : lm.asMap().entrySet()) {
			if(max < ent.getValue().size()){
				max = ent.getValue().size();
			}
			if(!first){
				bw.write("\t");
				
			}
			first = false;
			bw.write(prefix);
			bw.write(ent.getKey().toString());
		}
		bw.write("\n");
		for(int i = 0; i < max; i++){
			first = true;
			for(A a : lm.keySet()){
				List<B> b = lm.get(a);
				if(!first){
					bw.write("\t");
				}
				first = false;
				if(b.size() < i){
					bw.write(b.get(i).toString());
				}
				else{
					bw.write(na);
				}
			}
			bw.write("\n");
		}
		bw.flush();
		bw.close();
	}

	public static void saveToFile(Map<?, ?> data, String file, String delimiter)
			throws Exception {
		DataOutputStream out = new DataOutputStream(new FileOutputStream(file));
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out));
		for (Map.Entry<?, ?> ent : data.entrySet()) {
			bw.write(ent.getKey().toString());
			bw.write(delimiter);
			Object object = ent.getValue();
			if ((object instanceof Collection)) {
				writeRow(bw, (Collection) object, delimiter);
			} else {
				bw.write(object.toString());
			}
			bw.write("\n");
		}
		bw.flush();
		out.close();
	}
	
	public static void  saveToFileEx(Map<String, Collection<String>> map, String file, String delimiter) throws Exception {
		DataOutputStream out = new DataOutputStream(new FileOutputStream(file));
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out));
		for (Map.Entry<String, Collection<String>> ent : map.entrySet()) {
			for(String o : ent.getValue()){
				bw.write(ent.getKey().toString());
				bw.write(delimiter);
				bw.write(o.toString());
				bw.write("\n");
			}
		}
		bw.flush();
		out.close();
	}
	
	public static void  saveToFileEx(Map<String, Collection<String>> map, String file, String delimiter, Map<String, String> append) throws Exception {
		DataOutputStream out = new DataOutputStream(new FileOutputStream(file));
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(out));
		for (Map.Entry<String, Collection<String>> ent : map.entrySet()) {
			for(String o : ent.getValue()){
				bw.write(ent.getKey().toString());
				bw.write(delimiter);
				bw.write(o.toString());
				if(append != null){
					bw.write(delimiter);
					String app = append.get(o);
					if(app != null){
						bw.write(app);
					}
				}
				bw.write("\n");
			}
		}
		bw.flush();
		out.close();
	}

	public static void saveToFile(Collection<?> data, String file)
			throws Exception {
		saveToFile(data, file, "\t");
	}
	
	public static void saveToFile(StringBuffer sb, String file)throws Exception {
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		bw.write(sb.toString());
		bw.flush();
		bw.close();
	}
	
	public static void saveToFile(StringBuffer sb, String file, boolean append)throws Exception {
		BufferedWriter bw = new BufferedWriter(new FileWriter(file, append));
		bw.write(sb.toString());
		bw.flush();
		bw.close();
	}
	
	public static Set<ArrayKey<String>> readEdgeList(String file, int col1, int col2, String delimiter) throws Exception{
		Set<ArrayKey<String>> set = new HashSet<ArrayKey<String>>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			String[] data = strLine.split(delimiter);
			String [] arr = new String[]{data[col1], data[col2]};
			Arrays.sort(arr);
			set.add(new ArrayKey<String>(arr));
		}
		br.close();
		return set;
	}
	
	public static ListMultimap<ArrayKey<String>, String> readEdges(String file, int col1, int col2, String delimiter) throws Exception{
		ListMultimap<ArrayKey<String>, String> map = ArrayListMultimap.create();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			String[] data = strLine.split(delimiter);
			String [] arr = new String[]{data[col1], data[col2]};
			Arrays.sort(arr);
			ArrayKey<String> key = new ArrayKey<String>(arr);
			for(int i = 0; i < data.length; i++){
				if(i != col1 && i != col2){
					map.put(key, data[i]);	
				}
			}
		}
		br.close();
		return map;
	}
	
	public static ListMultimap<ArrayKey<String>, String> combine(Collection<ArrayKey<String>> n1, Collection<ArrayKey<String>> n2){
		ListMultimap<ArrayKey<String>, String> map = ArrayListMultimap.create();
		Set<ArrayKey<String>> set = new HashSet<ArrayKey<String>>();
		set.addAll(n1);
		set.addAll(n2);
		for(ArrayKey<String> key : set){
			map.put(key, n1.contains(key)?"1":"0");
			map.put(key, n2.contains(key)?"1":"0");	
		}
		return map;
	}
	
	public static void append(ListMultimap<ArrayKey<String>, String> map, Collection<ArrayKey<String>> n2){
		Set<ArrayKey<String>> set = new HashSet<ArrayKey<String>>();
		set.addAll(map.keySet());
		set.addAll(n2);
		for(ArrayKey<String> key : set){
			map.put(key, map.containsKey(key)?"1":"0");
			map.put(key, n2.contains(key)?"1":"0");	
		}
	}

	public static void saveToFile(Map<?, ?> data, String file) throws Exception {
		saveToFile(data, file, "\t");
	}

	public static ListMultimap<String, String> readMap(String file, String delimiter) throws Exception {
		return readMap(file, delimiter, ArrayListMultimap.<String, String>create());
	}

	public static ListMultimap<String, String> readMap(String file, String delimiter, ListMultimap<String, String> map) throws Exception {
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			String[] data = strLine.split(delimiter);
			if (data.length > 1) {
				for (int i = 1; i < data.length; i++) {
					map.put(data[0].trim(), data[i]);
				}
			}
		}
		in.close();
		return map;
	}
	
	public static ListMultimap<Integer, String> readMapInt(String file, String delimiter) throws Exception {
		return readMapInt(file, delimiter, ArrayListMultimap.<Integer, String>create());
	}

	public static ListMultimap<Integer, String> readMapInt(String file, String delimiter, ListMultimap<Integer, String> map) throws Exception {
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			String[] data = strLine.split(delimiter);
			if (data.length > 1) {
				for (int i = 1; i < data.length; i++) {
					map.put(Integer.valueOf(data[0].trim()), data[i]);
				}
			}
		}
		in.close();
		return map;
	}
	
	public static <A, B> void  saveListMultimap(ListMultimap<A, B> map, String file, String delimiter) throws Exception{
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		for (Entry<A, Collection<B>> ent : map.asMap().entrySet()) {
			bw.write(ent.getKey().toString());
			for(B b : ent.getValue()){
				bw.write(delimiter);
				bw.write(b.toString());
			}
			bw.write("\n");
		}
		bw.flush();
		bw.close();
	}
	
	public static ListMultimap<Double, String> readMapDbl(String file, String delimiter) throws Exception {
		return readMapDbl(file, delimiter, ArrayListMultimap.<Double, String>create());
	}

	public static ListMultimap<Double, String> readMapDbl(String file, String delimiter, ListMultimap<Double, String> map) throws Exception {
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			String[] data = strLine.split(delimiter);
			if (data.length > 1) {
				for (int i = 1; i < data.length; i++) {
					map.put(Double.valueOf(data[0].trim()), data[i]);
				}
			}
		}
		in.close();
		return map;
	}
	

	
	public static Map<String, Double> readSimpleDblMap(String file, String delimiter, int col1, int col2) throws Exception {
		Map<String, Double> map = new HashMap<String, Double>();
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			String[] data = strLine.split(delimiter);
			if (data.length > 1) {
				String a = data[col1].trim();
				String b = data[col2].trim();
				if(!a.equals("") && !b.equals("")){
					try{
						map.put(a, Double.valueOf(b.trim()));
					}
					catch(NumberFormatException e){
						map.put(a, Double.NaN);
					}
				}
			}
		}
		in.close();
		return map;
	}

	public static Map<String, String> readSimpleMap(String file, String delimiter, int col1, int col2) throws Exception {
		Map<String, String> map = new HashMap<String, String>();
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			String[] data = strLine.split(delimiter);
			if (data.length > 1) {
				String a = data[col1].trim();
				String b = data[col2].trim();
				if(!a.equals("") && !b.equals("")){
					map.put(a, b);
				}
				
			}
		}
		in.close();
		return map;
	}

	public static Map<String, String> readSimpleMap1(String file, String delimiter, int col1, int col2) throws Exception {
		Map<String, String> map = new HashMap<String, String>();
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			String[] data = strLine.split(delimiter);
			if (data.length > 1) {
				String a = data[col1].trim();
				String b = data[col2].trim();
				if(!a.equals("") && !b.equals("")){
					map.put(a, b);
				}
				
			}
		}
		in.close();
		return map;
	}
	
	public static Map<String, String> readSimpleMap(String file, String delimiter) throws Exception {
		return readSimpleMap(file, delimiter, 0, 1);
	}

	public static List<List<String>> readList(String file, String delimiter) throws Exception {
		List<List<String>> list = new LinkedList<List<String>>();
		DataInputStream in = new DataInputStream(new FileInputStream(file));
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String strLine;
		while ((strLine = br.readLine()) != null) {
			list.add(Arrays.asList(strLine.split(delimiter)));
		}
		in.close();
		return list;
	}
	
	public static class AnnotationComparator implements Comparator<String[]>{
		private final int [] pos;
		
		public AnnotationComparator(int ... pos){
			this.pos = pos;
		}
		
		int v = 0;
		public int compare(String[] o1, String[] o2) {
			for(int i : this.pos){
				int v = o1[i].compareTo(o2[i]);
				if(v != 0){
					return v;
				}
			}
			return v;
		}
	}
	
	public static class AnnotationComparatorInt implements Comparator<String[]>{
		private final int [] pos;
		
		public AnnotationComparatorInt(int ... pos){
			this.pos = pos;
		}
		
		int v = 0;
		public int compare(String[] o1, String[] o2) {
			for(int i : this.pos){
				int v = Integer.valueOf(o1[i]).compareTo(Integer.valueOf(o2[i]));
				if(v != 0){
					return v;
				}
			}
			return v;
		}
	}
}
