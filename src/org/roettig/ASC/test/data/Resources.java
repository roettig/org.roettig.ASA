package org.roettig.ASC.test.data;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class Resources
{
	public static String readString(InputStream in) throws IOException
	{
		BufferedReader rd = new BufferedReader(new InputStreamReader(in));
		
		StringBuffer sb = new StringBuffer();
		
		String line="";
		while( (line=rd.readLine())!=null )
		{
			sb.append(line+ System.getProperty("line.separator"));
		}
		
		return sb.toString();
	}
}
