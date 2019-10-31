package uk.ac.rothamsted.phinets;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;

public class BiomartDownloader {
	private static final int CONNECTION_TIMEOUT = 3000;
	private static final int READ_TIMEOUT = 5000;
	private static final String sc = "scerevisiae";
	private static final String sp = "spombe";
	
	public static final String orthSc = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"fungi_mart\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.7\" ><Dataset name = \"&&&&_eg_gene\" interface = \"default\" ><Attribute name = \"ensembl_transcript_id\" /><Attribute name = \"####_eg_gene\" /><Attribute name = \"####_eg_homolog_is_tree_compliant\" /><Attribute name = \"####_eg_orthology_type\" /></Dataset></Query>";
	public static final String goUrl = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"fungi_mart\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.7\" ><Dataset name = \"&&&&_eg_gene\" interface = \"default\" ><Attribute name = \"ensembl_transcript_id\" /><Attribute name = \"go_accession\" /></Dataset></Query>";
	public static final String phiUrl = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"fungi_mart\" formatter = \"TSV\" header = \"0\" uniqueRows = \"0\" count = \"\" datasetConfigVersion = \"0.7\" ><Dataset name = \"&&&&_eg_gene\" interface = \"default\" ><Attribute name = \"ensembl_transcript_id\" /><Attribute name = \"phibase_phenotype\" /></Dataset></Query>";
	
	public static void main(String[] args) throws Exception {
		//generateConfig();
		download();
	}
	
	public static final void download() throws Exception{
		String speciesListFile = "D:\\projects\\path_net_combined\\domains\\species_list.tab";
		String orthoDir = "D:\\projects\\path_net_combined\\orthologs\\";
		String goDir = "D:\\projects\\path_net_combined\\goa\\";
		String phiDir = "D:\\projects\\path_net_combined\\phibase\\";
		BufferedReader br = new BufferedReader(new FileReader(speciesListFile));
		String line;
		
		/*
		*/
		while ((line = br.readLine()) != null) { 
			String[] data = line.split("\t");
			for(String orth : new String[]{sp, sc}){
				try{
					downloadOrthologs(data[0]+"_vs_"+orth+".tab", orthoDir, orth, data[2]);
				}
				catch(RuntimeException e){
					e.printStackTrace();
				}
				catch(Exception e){
					e.printStackTrace();
				}
			}
		}
		br.close();
		
		br = new BufferedReader(new FileReader(speciesListFile));
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
				try{
					downloadGoa(data[0]+"_goa.tab", goDir, data[2]);
				}
				catch(RuntimeException e){
					e.printStackTrace();
				}
				catch(Exception e){
					e.printStackTrace();
				}
			
		}
		br.close();
		
		br = new BufferedReader(new FileReader(speciesListFile));
		while ((line = br.readLine()) != null) {
			String[] data = line.split("\t");
				try{
					downloadPhibase(data[0]+"_phi.tab", phiDir, data[2]);
				}
				catch(RuntimeException e){
					e.printStackTrace();
				}
				catch(Exception e){
					e.printStackTrace();
				}
			
		}
		br.close();
		
	}
	
	public static void generateConfig() throws Exception{
		String dir = "D:\\projects\\path_net_combined\\domains\\";
		BufferedWriter bw = new BufferedWriter(new FileWriter(dir+"species_list.tab"));
		for(File f : new File(dir).listFiles()){
			System.err.println(f.getName());
			if(!f.getName().startsWith("0")){
				continue;
			}
			String [] tmp = f.getName().split("\\.")[0].split("_");
			bw.write(tmp[1]+"_"+tmp[2]);
			bw.write("\t");
			bw.write(f.getName());
			bw.write("\t");
			bw.write(tmp[1].substring(0, 1).toLowerCase()+tmp[2]);
			bw.write("\n");
			bw.flush();

		}
		bw.close();

	}
	
	public static String downloadOrthologs(String file, String dir, String ortho, String species) throws Exception {
	   //String species = "fgraminearum";
	    String link = orthSc.replace("&&&&", species);
	    link = link.replaceAll("####", ortho);
	    URL url = new URL("http://fungi.ensembl.org/biomart/martservice?query="+URLEncoder.encode(link, "UTF-8"));
	    final URLConnection connection = url.openConnection();
	    connection.setReadTimeout(READ_TIMEOUT);
		connection.setConnectTimeout(CONNECTION_TIMEOUT);

		InputStream is = connection.getInputStream();
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(is));
		String line = reader.readLine();
		String fileName = dir+file;
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		
		while ((line = reader.readLine()) != null) {
			String[] data = line.split("\t");
			if(data.length < 2 || data[0].trim().isEmpty() || data[1].trim().isEmpty()){
				continue;
			}
			bw.write(line+"\n");
		}
		bw.flush();
		bw.close();
		
		return fileName;
	}
	
	public static String downloadGoa(String file, String dir, String species) throws Exception {
		   //String species = "fgraminearum";
		    String link = goUrl.replace("&&&&", species);
		    URL url = new URL("http://fungi.ensembl.org/biomart/martservice?query="+URLEncoder.encode(link, "UTF-8"));
		    final URLConnection connection = url.openConnection();
		    connection.setReadTimeout(READ_TIMEOUT);
			connection.setConnectTimeout(CONNECTION_TIMEOUT);

			InputStream is = connection.getInputStream();
			
			BufferedReader reader = new BufferedReader(new InputStreamReader(is));
			String line = reader.readLine();
			String fileName = dir+file;
			BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
			
			while ((line = reader.readLine()) != null) {
				String[] data = line.split("\t");
				if(data.length < 2 || data[0].trim().isEmpty() || data[1].trim().isEmpty()){
					continue;
				}
				bw.write(line+"\n");
			}
			bw.flush();
			bw.close();
			
			return fileName;
		}
	
	public static String downloadPhibase(String file, String dir, String species) throws Exception {
		   //String species = "fgraminearum";
		    String link = phiUrl.replace("&&&&", species);
		    URL url = new URL("http://fungi.ensembl.org/biomart/martservice?query="+URLEncoder.encode(link, "UTF-8"));
		    final URLConnection connection = url.openConnection();
		    connection.setReadTimeout(READ_TIMEOUT);
			connection.setConnectTimeout(CONNECTION_TIMEOUT);

			InputStream is = connection.getInputStream();
			
			BufferedReader reader = new BufferedReader(new InputStreamReader(is));
			String line = reader.readLine();
			String fileName = dir+file;
			BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
			
			while ((line = reader.readLine()) != null) {
				String[] data = line.split("\t");
				if(data.length < 2 || data[0].trim().isEmpty() || data[1].trim().isEmpty()){
					continue;
				}
				bw.write(line+"\n");
			}
			bw.flush();
			bw.close();
			
			return fileName;
		}
}
