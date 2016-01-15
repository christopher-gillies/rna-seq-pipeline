package org.kidneyomics.rnaseq;

import java.io.IOException;
import java.util.List;

public interface BAMInfoService {
	List<BAMInfo> getBAMInfosFromBamList(String bamList) throws IOException;
}
