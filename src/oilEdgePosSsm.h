#ifndef __EDGEPOSSSM_H__
#define __EDGEPOSSSM_H__

#define SNAME_OILEDGEPOS "oilEdgePos"

typedef struct{
	int dist;      // 溶接エッジまでの距離[mm]
	double theta;  // 溶接エッジの角度[rad]
	int step;      // サイド溶接段差検知結果( 0:なし, 1:左のみ, 2:右のみ, 3:両サイド )
	int status;    // 検知結果のステータス (1:成功, 0:失敗)
}OilEdgePos;

#endif

