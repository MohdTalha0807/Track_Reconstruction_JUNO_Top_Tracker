#ifndef CHANNEL_ID_TOOLS_H
#define CHANNEL_ID_TOOLS_H

unsigned int getDetectorID(const unsigned int ch_id){
	return (ch_id>>24) & 0xFF;
}

unsigned int getLayerID(const unsigned int ch_id){
	return (ch_id>>16) &  0x7;
}

unsigned int getColumnID(const unsigned int ch_id){
	return (ch_id>>13) &  0x7;
}

unsigned int getRowID(const unsigned int ch_id){
	return (ch_id>>10) &  0x7;
}

unsigned int getWallID(const unsigned int ch_id){
	return (ch_id>>10) &  0x1FF;
}

unsigned int getWallHID(const unsigned int ch_id){
	return getLayerID(ch_id)*100 + getColumnID(ch_id)*10 + getRowID(ch_id);
}

unsigned int getPMTID(const unsigned int ch_id){
	return (ch_id>> 6) &  0xF;
}

unsigned int getStripID(const unsigned int ch_id){
	return (ch_id    ) & 0x3F;
}

void decodeChannelID(const unsigned int ch_id, unsigned int & detector, unsigned int & layer, unsigned int & column, unsigned int & row, unsigned int & PMT, unsigned int & strip){
	detector = (ch_id>>24) & 0xFF;
	layer    = (ch_id>>16) &  0x7;
	column   = (ch_id>>13) &  0x7;
	row      = (ch_id>>10) &  0x7;
	PMT      = (ch_id>> 6) &  0xF;
	strip    = (ch_id    ) & 0x3F;
}

unsigned int makeChannelID(unsigned int detector, unsigned int layer, unsigned int column, unsigned int row, unsigned int PMT, unsigned int strip){
	unsigned int ch_id=0;
	ch_id += (detector & 0xFF)<<24;
	ch_id += (layer    &  0x7)<<16;
	ch_id += (column   &  0x7)<<13;
	ch_id += (row      &  0x7)<<10;
	ch_id += (PMT      &  0xF)<<6;
	ch_id += (strip    & 0x3F)    ;
	return ch_id;
}

bool isInTT(const unsigned int channelID){
	unsigned int detector, layer,  column, row, PMT, strip;
	decodeChannelID(channelID, detector, layer,  column, row, PMT, strip);
	if(detector!=0x30) return false;
	if(layer>=3) return false;
	if(column>=3) return false;
	if(row>=7) return false;
	if(PMT>=16) return false;
	if(strip>=64) return false;
	return true;
}

bool isXPlane(const unsigned int channelID){
	unsigned int PMT = getPMTID(channelID) ;
	switch(PMT){
		case 0:
		case 1:
		case 2:
		case 3:
		case 8:
		case 9:
		case 10:
		case 11:
			return true;
		default:
			return false;
	}
}

#endif

