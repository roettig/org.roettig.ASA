package org.roettig.ASC.asc.interceptors;

import org.roettig.ASC.asc.ASAFlow;
import org.roettig.ASC.asc.ASAFlow.ASAFlowContext;

public interface ASAFlowInterceptor
{
	public void intercept(ASAFlow flow, ASAFlowContext ctx);
}
