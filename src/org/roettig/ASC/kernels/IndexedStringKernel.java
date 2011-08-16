package org.roettig.ASC.kernels;

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import org.roettig.MLToolbox.base.instance.DualInstance;
import org.roettig.MLToolbox.base.parameter.Parameter;
import org.roettig.MLToolbox.kernels.KernelFunction;
import org.roettig.MLToolbox.kernels.RBFKernel;

public class IndexedStringKernel extends KernelFunction<DualInstance<String>>
{
	private List<Integer> fixedIdx = null;
	private KernelFunction<DualInstance<String>> k_fun;
		
    public IndexedStringKernel(KernelFunction<DualInstance<String>> k_fun)
    {
        fixedIdx = new ArrayList<Integer>();
        this.k_fun = k_fun;
    }
    
    public IndexedStringKernel(KernelFunction<DualInstance<String>> k_fun, int idx)
    {
    	this(k_fun);
        addFixedIndex(idx);
    }

    public void addFixedIndex(int i)
    {
        fixedIdx.add( i );
    }
    
    public void clearIndices()
    {
    	fixedIdx.clear();
    }
    
	@Override
	public double compute(DualInstance<String> x, DualInstance<String> y) throws Exception
	{
        String xs = x.getPayload(); 
        String ys = y.getPayload();
     
        if(xs.length()!=ys.length())
        	throw new Exception("error evaluating kernel function "+this.getClass().getCanonicalName()+" : supplied strings have different lengths");
        
        double s = 0.0;
        
        for(int i=0;i<xs.length();i++)
        {
           if(fixedIdx.contains(i))
           {
        	   s += k_fun.compute( new DualInstance<String>( null, xs.substring(i,i+1)), new DualInstance<String>( null, ys.substring(i,i+1)));
        	   //System.out.print(s+" ");
           }
        }
        
        return s;
	}
}