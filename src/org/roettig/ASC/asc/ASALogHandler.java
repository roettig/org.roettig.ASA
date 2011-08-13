/**
 * 
 */
package org.roettig.ASC.asc;

import java.util.logging.Handler;
import java.util.logging.LogRecord;

/**
 * @author roettig
 *
 */
public class ASALogHandler extends Handler
{
	private ASAFlow flow;

	public ASALogHandler(ASAFlow _flow)
	{
		flow = _flow;
	}

	@Override
	public void close() throws SecurityException
	{
	}

	@Override
	public void flush()
	{
	}

	@Override
	public void publish(LogRecord arg0)
	{	
		flow.log("",arg0.getMessage().toString());
	}
}
