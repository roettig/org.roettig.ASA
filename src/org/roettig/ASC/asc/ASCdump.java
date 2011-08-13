/**
 * 
 */
package org.roettig.ASC.asc;

import java.io.FileNotFoundException;
import java.io.IOException;

/**
 * @author roettig
 *
 */
public class ASCdump
{

	/**
	 * @param args
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException, IOException
	{
		String jobid = "/tmp/job.ini";//args[0];
		ASAJob jd = null;
		boolean failed = false;
		try
		{
			jd= ASAJob.load(jobid);
		}
		catch(Exception e)
		{
			failed = true;
		}
		if(failed)
		{
			jd = JobFileGenerator.read(jobid);    
		}
		System.out.println("alnmethod="+jd.getMSAmethod());
	}

}
