#ifndef MANTID_ISISREFLECTOMETRY_QWIDGETGROUP_H
#define MANTID_ISISREFLECTOMETRY_QWIDGETGROUP_H
#include <cstddef>
#include <array>
#include <QWidget>
namespace MantidQt {
namespace CustomInterfaces {
template <std::size_t N> class QWidgetGroup {
public:
  QWidgetGroup() : m_widgets() {}
  explicit QWidgetGroup(std::array<QWidget *, N> const &widgets)
      : m_widgets(widgets) {}

  void enable() {
    for (auto *widget : m_widgets)
      widget->setEnabled(true);
  }

  void disable() {
    for (auto *widget : m_widgets)
      widget->setEnabled(false);
  }

private:
  std::array<QWidget *, N> m_widgets;
};

template <typename... Ts>
QWidgetGroup<sizeof...(Ts)> makeQWidgetGroup(Ts... widgets) {
  return QWidgetGroup<sizeof...(Ts)>(
      std::array<QWidget *, sizeof...(Ts)>({widgets...}));
}
}
}
#endif // MANTID_ISISREFLECTOMETRY_QWIDGETGROUP_H
