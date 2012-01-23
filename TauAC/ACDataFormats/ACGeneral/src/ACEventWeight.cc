#include "../interface/ACEventWeight.h"

ClassImp(ACEventWeight)

ACEventWeight::ACEventWeight(double weight) {
    weight_ = weight;
}
ACEventWeight::ACEventWeight() {
    weight_ = -1.0;
}
