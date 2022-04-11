/**
 *
 *
 * @file ContainerContainer.h
 * @brief 
 * @author clonker
 * @date 4/8/22
 */
#pragma once

namespace ctiprd::cpu {

template<typename... Containers>
class ContainerContainer {
public:

    auto begin() {

    }

    auto end() {

    }

    auto size() {

    }

private:
    std::tuple<Containers...> containers {};
};

}
