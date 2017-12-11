#pragma once

enum class Species : int { mReactant, mProduct, mCation, mPotential, sReactant, sProduct, sCation, sPotential, sAnion, Count };

enum class boundaries {mbottom, mbulk, mtop, sbottom, sbulk, stop};

enum class SignalStateSwitch {endOfPositive, endOfNegative, others};