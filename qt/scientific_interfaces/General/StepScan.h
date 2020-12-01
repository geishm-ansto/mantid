// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#pragma once

//----------------------
// Includes
//----------------------
#include "MantidAPI/IAlgorithm.h"
#include "MantidAPI/MatrixWorkspace_fwd.h"
#include "MantidQtWidgets/Common/AlgorithmRunner.h"
#include "MantidQtWidgets/Common/UserSubWindow.h"
#include "ui_StepScan.h"

namespace MantidQt {
namespace CustomInterfaces {

class StepScan : public API::UserSubWindow {
  Q_OBJECT

public:
  /// The name of the interface as registered into the factory
  static std::string name() { return "Step Scan Analysis"; }
  // This interface's categories.
  static QString categoryInfo() { return "General"; }

  explicit StepScan(QWidget *parent = nullptr);
  ~StepScan() override;

signals:
  void logsAvailable(const Mantid::API::MatrixWorkspace_const_sptr & /*_t1*/);
  void logsUpdated(const Mantid::API::MatrixWorkspace_const_sptr & /*_t1*/);
  void updatePlot(const QString & /*_t1*/);

private slots:
  void triggerLiveListener(bool checked);
  void startLiveListenerComplete(bool error);
  void loadFile(bool async = true);
  void loadFileComplete(bool error);
  void launchInstrumentWindow();
  void fillPlotVarCombobox(const Mantid::API::MatrixWorkspace_const_sptr &ws);
  void expandPlotVarCombobox(const Mantid::API::MatrixWorkspace_const_sptr &ws);
  void fillNormalizationCombobox();
  void runStepScanAlg();
  bool runStepScanAlgLive(const std::string &stepScanProperties);

  void updateForNormalizationChange();
  void generateCurve(const QString &var);

  void helpClicked();

private:
  void initLayout() override;
  void startLiveListener();
  bool mergeRuns();
  void setupOptionControls();
  void clearNormalizationCombobox();
  Mantid::API::IAlgorithm_sptr setupStepScanAlg();
  void cleanupWorkspaces();
  void plotCurve();

  void handleAddEvent(Mantid::API::WorkspaceAddNotification_ptr pNf);
  void handleReplEvent(Mantid::API::WorkspaceAfterReplaceNotification_ptr pNf);
  void addReplaceObserverOnce();
  void checkForMaskWorkspace(const std::string &wsName);
  void checkForResultTableUpdate(const std::string &wsName);
  void checkForVaryingLogs(const std::string &wsName);

  Ui::StepScan m_uiForm; ///< The form generated by Qt Designer
  std::string m_inputWSName, m_tableWSName, m_plotWSName;
  QString m_inputFilename;
  const std::string m_instrument; ///< The default instrument (for live data)

  API::AlgorithmRunner
      *m_algRunner; ///< Object for running algorithms asynchronously
  Poco::NObserver<StepScan, Mantid::API::WorkspaceAddNotification>
      m_addObserver;
  Poco::NObserver<StepScan, Mantid::API::WorkspaceAfterReplaceNotification>
      m_replObserver;
#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
  std::optional<int> m_fignum;
#endif
  bool m_replaceObserverAdded;
};

} // namespace CustomInterfaces
} // namespace MantidQt
