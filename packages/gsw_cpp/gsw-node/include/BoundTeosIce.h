//
// Created by blake on 9/26/24.
//

#ifndef BOUNDTEOSICE_H
#define BOUNDTEOSICE_H

#include <napi.h>

#include <TeosIce.h>

class BoundTeosIce : public Napi::ObjectWrap<BoundTeosIce>
{
public:
    static Napi::Object Init(Napi::Env env, Napi::Object exports);

    explicit BoundTeosIce(const Napi::CallbackInfo& info);

    Napi::Value gsw_cp_ice(const Napi::CallbackInfo& info);

private:
    TeosIce m_ice;
};



#endif //BOUNDTEOSICE_H
