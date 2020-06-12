#include "MarketParam.h"

double MarketParam::get_Lvol(int i_t, int i_k) const
{
	return vol.get_Lvol(i_t, i_k);
}
