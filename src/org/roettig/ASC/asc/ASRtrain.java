package org.roettig.ASC.asc;

import java.io.File;

import org.roettig.ASC.asc.interceptors.AlignmentInterceptor;
import org.roettig.SequenceTools.MSA;


/**
 * @author roettig
 *
 */
public class ASRtrain
{

	public static void main(String[] args) throws Exception
	{
		ASAJob jd = null;
		boolean failed = false;
		try
		{
			jd= ASAJob.load(args[0]);
		}
		catch(Exception e)
		{
			System.out.println("Jobfile "+args[0]+" seems to be no binary jobfile. Trying text-mode.");
			failed = true;
		}
		if(failed)
		{
			jd = JobFileGenerator.read(args[0]);    
		}

		jd.setOutputDir(".");

		ASAFlow ascflow = new ASR(jd);


		File f = new File(jd.getOutputDir()+"/core.afa");
		if(f.exists())
		{
			System.out.println("Restoring core msa : "+f.getAbsolutePath());
			MSA msa = MSA.loadFromFile(f.getAbsolutePath());
			ascflow.setAlignmentInterceptor( new AlignmentInterceptor(new MockAlignmentBuilder(msa)) );
		}
		else
		{
			if(jd.hasThreeDCoffeeId())
			{
				System.out.println("Using stored 3DCOFFEE core msa");
				ascflow.setAlignmentInterceptor( new AlignmentInterceptor(new ThreeDCoffeeAlignmentBuilder(jd.getThreeDCoffeeId())) );
			}
			else
			{
				System.out.print("Building new core msa using ["+jd.getMSAmethod()+"]");
				if(jd.getMSAmethod().equalsIgnoreCase("3dcoffee"))
				{
					ascflow.setAlignmentInterceptor( new AlignmentInterceptor(new ThreeDCoffeeAlignmentBuilder()) );
					System.out.println(" 3dcoffee");
				}
				else
				{
					ascflow.setAlignmentInterceptor( new AlignmentInterceptor(new MuscleAlignmentBuilder()) );
					System.out.println(" muscle");
				}
			}
		}

		ascflow.run();

		ASAModel mdl = new ASAModel(ascflow.getModel(),ascflow.getMSA());
		mdl.setASCresidues(jd.getASCresidues());
		mdl.save("ASRmodel");
		jd.save("ASAdata");
	}

}
