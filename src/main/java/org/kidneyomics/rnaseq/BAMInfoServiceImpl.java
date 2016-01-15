package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.springframework.stereotype.Component;

@Component(value="bamInfoService")
public class BAMInfoServiceImpl implements BAMInfoService {

	@Override
	public List<BAMInfo> getBAMInfosFromBamList(String bamListPath) throws IOException {
		
		File bamList = new File(bamListPath);
		
		if(!bamList.exists()) {
			throw new IllegalArgumentException("bamList "+ bamList.getAbsolutePath() + " not found!");
		}
		
		List<String> lines = FileUtils.readLines(bamList);
		List<BAMInfo> infos = BAMInfo.getBAMInfoFromLines(lines);
		
		return infos;
	}

}
