/**
 * 
 */
package org.roettig.ASC.asc;

import org.roettig.SequenceTools.MSA;
import org.roettig.SequenceTools.SequenceSet;

/**
 * @author roettig
 *
 */
public interface AlignmentBuilder
{
    public MSA build(SequenceSet seqs) throws Exception;
}
