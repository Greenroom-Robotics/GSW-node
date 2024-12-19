//
// Created by blake on 9/26/24.
//

#ifndef BOUNDTEOSSEA_H
#define BOUNDTEOSSEA_H


#include <napi.h>

#include <TeosSea.h>


class BoundTeosSea : public Napi::ObjectWrap<BoundTeosSea>
{
public:
    static Napi::Object Init(Napi::Env env, Napi::Object exports);

    explicit BoundTeosSea(const Napi::CallbackInfo& info);

    Napi::Value gsw_c_from_sp(const Napi::CallbackInfo& info);

    Napi::Value gsw_sound_speed(const Napi::CallbackInfo& info);

private:
    TeosSea m_sea;
};


#endif //BOUNDTEOSSEA_H
