/**
 * 
 */
package org.roettig.ASC.asc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Vector;

import org.roettig.PDBTools.ResidueLocator;

/**
 * @author roettig
 *
 */
public class JobFileGenerator
{
	private static List<String> seqfiles   = new Vector<String>();
	private static List<String> classnames = new Vector<String>();
	private static List<ResidueLocator> reslocs = new Vector<ResidueLocator>();
	private static Map<String, List<Double>> params = new HashMap<String,List<Double> >();

	public static ASAJob read(String filename) throws FileNotFoundException, IOException
	{
		ASAJob ret = new ASAJob();
		Properties props = new Properties();
		props.load(new FileReader(new File(filename)) );

		Map<Integer,String> files = new HashMap<Integer,String>();
		Map<Integer,String> names = new HashMap<Integer,String>();

		for(Object keyO : props.keySet())
		{
			String key   = keyO.toString();
			String value = props.getProperty(key.toString());

			if(key.startsWith("seqfile"))
			{
				Integer idx = Integer.parseInt(key.replace("seqfile",""));
				files.put(idx, value);
				//seqfiles.add(value);
			}

			if(key.startsWith("classname"))
			{
				Integer idx = Integer.parseInt(key.replace("classname",""));
				names.put(idx, value);		
				//classnames.add(value);
			}

			if(key.toLowerCase().startsWith("pdbfile"))
			{
				ret.setPdbfile(value);
			}

			if(key.toLowerCase().startsWith("3dcoffeeid"))
			{
				ret.setThreeDCoffeeId(value);
			}

			if(key.toLowerCase().startsWith("outputdir"))
			{
				ret.setOutputDir(value);
			}

			if(key.toLowerCase().startsWith("distance"))
			{
				ret.setDistance(Double.parseDouble(value));
			}

			if(key.toLowerCase().startsWith("ninnerfolds"))
			{
				ret.setNInnerFolds(Integer.parseInt(value));
			}

			if(key.toLowerCase().startsWith("nouterfolds"))
			{
				ret.setNOuterFolds(Integer.parseInt(value));
			}

			if(key.toLowerCase().startsWith("resloc"))
			{
				reslocs.add( ResidueLocator.fromString(value) );
			}


			if(key.toLowerCase().startsWith("modeltype"))
			{
				ret.setModelType(value.toLowerCase());
			}

			if(key.toLowerCase().startsWith("kerneltype"))
			{
				ret.setKernelType(value.toLowerCase());
			}

			if(key.toLowerCase().startsWith("hyperparams"))
			{
				String toks[] = value.split(";");
				for(String tok: toks)
				{
					String toks2[] = tok.split(":");
					System.out.println("param "+toks2[0]);
					params.put(toks2[0], new Vector<Double>());
					String toks3[] = toks2[1].split(",");
					for(String val: toks3)
					{
						System.out.println("val:"+val);
						params.get(toks2[0]).add(Double.parseDouble(val));
					}
				}
			}

		}	

		ret.setHyperParams(params);
		ret.setReslocs(reslocs);

		for(Integer idx: files.keySet())
		{
			seqfiles.add(files.get(idx));
		}

		for(Integer idx: names.keySet())
		{
			classnames.add(names.get(idx));
		}


		ret.setFilenames(seqfiles);

		if(classnames.size()==0)
		{
			for(int i=0;i<seqfiles.size();i++)
			{
				classnames.add(String.format("class%d", i+1));
			}
		}

		ret.setClassnames(classnames);
		return ret;
	}

	public static void main(String[] args) throws FileNotFoundException, IOException
	{
		//String filename = args[0];
		//String filename = "/tmp/test.ini";
		//read(filename);
	}
}
